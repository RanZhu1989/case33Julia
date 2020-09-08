#NOTE  2020.9.6 v0.7 
#NOTE  Author：Ran Zhu @ School of Cyber Engineering and Science, SEU
#NOTE  此.jl函数可用于解决正常运行时配网运行优化问题和在一般故障时非孤岛运行条件下的自愈重构问题
#NOTE  通过设置startPoint 和 endPoint 可调节规划尺度，在此框架上可通过导入外部模型文件做成MPC规划器
#NOTE  **模型主要参考文献**：
#NOTE  1. 基于网络流模型的Brach Flow Model : 1989, M.E.Baran, TPS
#NOTE  2. Brach Flow Model 相角松弛法 ：2013, M.Farivar, TPS
#NOTE  3. Spanning Tree 拓扑约束 : 2012, J.A.Taylor, TPS
#NOTE  4. Single Commodity Flow 拓扑约束： 2013, R.A.Jabr, TPS
#NOTE  5. 开关操作次数约束： 2016, Mohammad, TPS
#NOTE  6. 判断节点和线路是否通电的约束驱动： 2019, Li, TSG
#NOTE  **计算机程序实现主要参考引用以下项目成果**：
#NOTE  1.MOSEK Project                  https://github.com/MOSEK/Mosek.jl
#NOTE  2.JuMP Modeling Language         https://arxiv.org/abs/1312.1431
#NOTE  3.Convex Programming Project     https://github.com/jump-dev/Convex.jl
#NOTE  4.MosekTools Project             https://github.com/JuliaOpt/MosekTools.jl
#NOTE  5.Julia Language Project @ MIT   https://github.com/JuliaLang/julia


                                               
#= __      _____  _ __| | __ (_)_ __    _ __  _ __ ___   __ _ _ __ ___ ___ ___ 
\ \ /\ / / _ \| '__| |/ / | | '_ \  | '_ \| '__/ _ \ / _` | '__/ _ / __/ __|
 \ V  V | (_) | |  |   <  | | | | | | |_) | | | (_) | (_| | | |  __\__ \__ \
  \_/\_/ \___/|_|  |_|\_\ |_|_| |_| | .__/|_|  \___/ \__, |_|  \___|___|___/
                                    |_|              |___/                  
 =#
 
import MathOptInterface
using DelimitedFiles
using Dates,TimesDates
using JuMP,MosekTools,Convex

#Set Current WorkPath
#ATTITION!  MUST make sure your workpath. 
cd(@__DIR__)
#Import Fuctions from .lib
include("fuctions_Lib.jl")

function dpSolverReconfiguraiton33Bus(typeOptimizer)
    #第一个参数：优化器选择  第二个参数
    currentStartTime=now()
    strCurrentStartTime=replace(replace(string(currentStartTime),"-"=>""),":"=>"")
    ###################建立JuMP模型##############################
    if typeOptimizer=="Mosek"
        bus33Reconfiguration=Model(Mosek.Optimizer)
    end
  
    ######################变量设置##################################
   
    #表示线路通断状态的alpha bin矩阵
    #ATTITION!这里从时间0开始，代表初始状态
    #包含故障线路
    @variable(bus33Reconfiguration,alpha[ijt in time0NodePair],Bin,base_name="线路通断状态alpha")

    #表示节点负载率情况
    @variable(bus33Reconfiguration,0<=loadRate[it in itPair]<=1,base_name="节点负载率loadRate")

    #用于构建辐射状拓扑约束的中SCF约束虚拟首端网络流变量apparentFictitiousFlow_{ijt}
    #不包含故障线路
    @variable(bus33Reconfiguration,apparentFictitiousFlow[ijt in ijtjit1Pair],base_name="虚拟网络流F")

    #用于构建开关次数约束的辅助变量lambda_{ijt}
    #ATTITION!这里从时间1开始
    #包含故障线路
    @variable(bus33Reconfiguration,lambda[ijt in time1NodePair],Bin,base_name="开关变位辅助变量lambda")

    #MT有功出力P_mt_{i t}
    @variable(bus33Reconfiguration,activePowerMT[it in tsMTalive],base_name="MT有功出力")
    
    #MT无功出力Q_mt_{i t}
    @variable(bus33Reconfiguration,reactivePowerMT[it in tsMTalive],base_name="MT无功出力")

    #注入有功\无功功率 P Q _{i t}
    @variable(bus33Reconfiguration,injectionActivePower[it in itPair],base_name="注入有功")
    @variable(bus33Reconfiguration,injectionReactivePower[it in itPair],base_name="注入无功")

    #首端有功\无功功率 P Q _{i j t}
    @variable(bus33Reconfiguration,apparentActivePower[ijt in time1NodePair],base_name="首端有功")
    @variable(bus33Reconfiguration,apparentReactivePower[ijt in time1NodePair],base_name="首端无功")

    #ATTITION! l_{i j t} 支路电流平方
    #每条支路最大载流量需要查lineData
    @variable(bus33Reconfiguration,lineSquCurrent[ijt in time1NodePair],base_name="线路电流幅值平方")

    #ATTITION!节点电压平方幅值 V_{i t}
    @variable(bus33Reconfiguration,nodeSquVoltage[it in itPair],base_name="节点电压幅值平方")

    #######################补充的变量取值范围#########################

    #故障线路在对应时刻一定是断开的
    for line in linesFault
        @constraint(bus33Reconfiguration,alpha[(line[1],line[2],line[3])]==0)  
    end
   
    #将初始状态的alpha(i,j,0)存入
    #ATTITION! 此初始状态指读取para_Line文件中init_State列的状态，每次solver的末尾必须更新init_State
    #TODO 一定记得在末尾加入更新init_State的语句
    for pair in lineInitStartPoint
        @constraint(bus33Reconfiguration,alpha[(pair[1],pair[2],0)]==pair[3]) 
    end

    #每条支路最大载流量需要查lineData
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,lineSquCurrent[(line[1],line[2],t)]<=alpha[(line[1],line[2],t)]*(maxCurrentLineDict[line])^2)
            @constraint(bus33Reconfiguration,lineSquCurrent[(line[1],line[2],t)]>=0)
        end
    end

    #MT出力约束
    for it in tsMTalive
        mtInf=mtInfDict[it[1]]
        @constraint(bus33Reconfiguration,mtInf[1]<=activePowerMT[it]<=mtInf[2])
        @constraint(bus33Reconfiguration,mtInf[3]<=reactivePowerMT[it]<=mtInf[4])
    end

    #节点电压约束
    for it in itPair
        @constraint(bus33Reconfiguration,nodeSquVoltage[it]<=highSquVNode)
        @constraint(bus33Reconfiguration,nodeSquVoltage[it]>=lowSquVNode)
    end

    #ATTITION! 这里假设为普通配电网形式，也可配置为微电网型，此时变电站出线节点电压设置为浮动
    #NOTE 如果是PV型DG 也可一并设置
    #将变电站出线节点电压设置为1.0pu~1.05pu
    for node in listSub
        for t in 1:points
            @constraint(bus33Reconfiguration,nodeSquVoltage[(node,t)]==voltageSquSub)
        end
    end

    #假设MT是具备黑启动能力的，即电压可以设置为1.0 pu.
    for it in tsMTalive
        @constraint(bus33Reconfiguration,nodeSquVoltage[it]==voltageBlackStartDG)
    end

    
    #注入功率：P_i=P发-P用
    #注意这里需要分以下情况：
    #1.纯负载节点或纯变电站节点;2.含PV的节点;3.含MT的节点
    #编译通过
    #NOTE 加入了削减负荷策略
    for i in rootFreeNodes
        for t in 1:points
            if (i,t) in tsMTalive
                @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==allocationDG[t,i][1]*pvActivePower[t,i]
                +allocationDG[t,i][2]*wtActivePower[t,i]+allocationDG[t,i][3]*activePowerMT[(i,t)]-consumeP[t,i]*loadRate[(i,t)])
                @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==allocationDG[t,i][1]*pvReactivePower[t,i]
                +allocationDG[t,i][2]*wtReactivePower[t,i]+allocationDG[t,i][3]*reactivePowerMT[(i,t)]-consumeQ[t,i]*loadRate[(i,t)])
            else
                @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==allocationDG[t,i][1]*pvActivePower[t,i]
                +allocationDG[t,i][2]*wtActivePower[t,i]-consumeP[t,i]*loadRate[(i,t)])
                @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==allocationDG[t,i][1]*pvReactivePower[t,i]
                +allocationDG[t,i][2]*wtReactivePower[t,i]-consumeQ[t,i]*loadRate[(i,t)])
            end
        end
    end
        

    #ATTITION! 针对传统配电网，变电站主变出线出力设置为大于等于0，作为MG时去掉此约束
    for i in listSub
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]>=0)
            @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]>=0)
        end
    end


    ########约束条件########################
    #网络流模型的功率平衡方程
    #编译通过
    for node in nodes
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(node,t)]
            ==sum(apparentActivePower[(j,k,t)] for (j,k) in findjkForwardNeighborPair(node,linesAlive,t) if j==node)
            -sum(apparentActivePower[(i,j,t)]-rLineDict[(i,j)]*lineSquCurrent[(i,j,t)] 
            for (i,j) in findijBackforwardNeighborPair(node,linesAlive,t) if j==node))
        
            @constraint(bus33Reconfiguration,injectionReactivePower[(node,t)]
            ==sum(apparentReactivePower[(j,k,t)] for (j,k) in findjkForwardNeighborPair(node,linesAlive,t) if j==node)
            -sum(apparentReactivePower[(i,j,t)]-xLineDict[(i,j)]*lineSquCurrent[(i,j,t)] 
            for (i,j) in findijBackforwardNeighborPair(node,linesAlive,t) if j==node))
        end
    end

    #节点电压联系方程
    #编译通过
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,nodeSquVoltage[(line[1],t)]-nodeSquVoltage[(line[2],t)]>=-bigM*(1-alpha[(line[1],line[2],t)])
            +2*(rLineDict[line]*apparentActivePower[(line[1],line[2],t)]+xLineDict[line]*apparentReactivePower[(line[1],line[2],t)])
            -((rLineDict[line])^2+(xLineDict[line]^2))*lineSquCurrent[(line[1],line[2],t)])

            @constraint(bus33Reconfiguration,nodeSquVoltage[(line[1],t)]-nodeSquVoltage[(line[2],t)]<=bigM*(1-alpha[(line[1],line[2],t)])
            +2*(rLineDict[line]*apparentActivePower[(line[1],line[2],t)]+xLineDict[line]*apparentReactivePower[(line[1],line[2],t)])
            -((rLineDict[line])^2+(xLineDict[line]^2))*lineSquCurrent[(line[1],line[2],t)])
        end
    end

    #首端功率P^2+Q^2=VI^2的旋转二阶锥约束
    #编译通过
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,[lineSquCurrent[(line[1],line[2],t)],0.5*nodeSquVoltage[(line[1],t)],
            apparentActivePower[(line[1],line[2],t)],apparentReactivePower[(line[1],line[2],t)]] in RotatedSecondOrderCone())
        end
    end

    #基于SCF生成树约束
    #Reference: R.A.Jabr,2013,TPS
    #sum_{i:i->j}F_ij,t +D_i,t =sum_{k:k->i}F_ki,t  i为非root节点
    #假设D_i,t==1
    #编译通过

    for node in rootFreeNodes
        for t in 1:points
            @constraint(bus33Reconfiguration,sum(apparentFictitiousFlow[(i,k,t)] for (i,k) in findjkForwardNeighborPair(node,linesAlive,t) if i==node)
            +1==sum(apparentFictitiousFlow[(k,i,t)] for (k,i) in findijBackforwardNeighborPair(node,linesAlive,t) if i==node))
        end
    end
    for line in linesAlive
        for t in 1:points
            @constraint(bus33Reconfiguration,(apparentFictitiousFlow[(line[1],line[2],t)])<=(bigM*alpha[(line[1],line[2],t)]))
        end
    end
    for line in linesAlive
        for t in 1:points
            @constraint(bus33Reconfiguration,(apparentFictitiousFlow[(line[1],line[2],t)])>=(-bigM*alpha[(line[1],line[2],t)]))
        end
    end
    for t in 1:points
        @constraint(bus33Reconfiguration,sum(alpha[(line[1],line[2],t)] for line in linesAlive)==numNode-1)
    end

  

    #开关操作次数约束
    #Reference: Mohammad 2016, TPS
    #需要加入去除 故障造成的上一步闭合的开关被断开的外部操作计数
    #编译通过
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,lambda[(line[1],line[2],t)]>=(alpha[(line[1],line[2],t)]-alpha[(line[1],line[2],t-1)]))
            @constraint(bus33Reconfiguration,lambda[(line[1],line[2],t)]>=(alpha[(line[1],line[2],t-1)]-alpha[(line[1],line[2],t)]))
        end
    end
    for t in 1:points
        @constraint(bus33Reconfiguration,(sum(lambda[(line[1],line[2],t)] for line in lines)
        -sum((alpha[(line[1],line[2],t-1)]-1)*1 for line in findTimeFalutLine(t,linesFault)))<=maxNS)  
    end

    ###################目标函数######################################

    #编译通过
    @objective(bus33Reconfiguration,Min,
    sum(mtInfDict[it[1]][5]*activePowerMT[(it[1],it[2])] for it in tsMTalive)
    +sum(sum(priceLoss*lineSquCurrent[(line[1],line[2],t)]*rLineDict[line] for line in lines)
    +sum(stGridPrice[t-1+startPoint]*injectionActivePower[(subNode,t)] for subNode in listSub)
    +sum(costSwitch*(sum(lambda[(line[1],line[2],t)] for line in lines)))
    +sum(failureLoss*(1-loadRate[(i,t)])*consumeP[t,i] for i in nodes) 
    for t in 1:points)
    -sum(sum((alpha[(line[1],line[2],t-1)]-1)*1*costSwitch for line in findTimeFalutLine(t,linesFault)) for t in 1:points if findTimeFalutLine(t,linesFault)!=())
    )
    
    #运行
    optimize!(bus33Reconfiguration)

    #Setting Path for Saving Report Files
    mkpath(".//report//")
    currentEndTime=now()
    strCurrentEndTime=replace(replace(string(currentEndTime),"-"=>""),":"=>"")
    strRunTime=parse(Float64,split(string(currentEndTime-currentStartTime)," ")[1])/1e4
    timeCell="_TimeCell from $startPoint to $endPoint "    
    reportPath=".//report//"
    resultPath= strCurrentStartTime * timeCell* "_简报.txt"
    resultSwitchOpearingPath=reportPath * timeCell* "_开关动作_" * strCurrentStartTime * resultPath
    resultPowerMTPath=reportPath *timeCell*"_MT出力_kW_" * resultPath
    resultNodeVoltagePath=reportPath * timeCell*"_节点电压_pu._" * strCurrentStartTime * resultPath
    resultLineCurrentPath=reportPath * timeCell*"_线路电流_A_" * strCurrentStartTime * resultPath
    resultLossPath=reportPath * timeCell*"_线路网损_kW_" * strCurrentStartTime * resultPath
    resultPowerSubstationPath=reportPath * timeCell*"_变电站主变出线出力_" * strCurrentStartTime *resultPath
    resultLoadEnergrized=reportPath * timeCell*"_负荷供电状态_" * strCurrentStartTime * resultPath

    #结果分析
    #NOTE 获取优化器结果标志位
    state=termination_status(bus33Reconfiguration)
    #NOTE 总报告
    fileResult=open(resultPath,"w")
    #NOTE 开关动作时间序列
    switchOperationFileResult=open(resultSwitchOpearingPath,"w")
    #NOTE MT出力调节时间序列
    powerMTFileResult=open(resultPowerMTPath,"w")
    #NOTE 节点电压幅值时间序列
    nodeVoltageFileResult=open(resultNodeVoltagePath,"w")
    #NOTE 线路电流时间序列
    lineCurrentFileResult=open(resultLineCurrentPath,"w")
    #NOTE 线路网损时间序列
    lineLossFileResult=open(resultLossPath,"w")
    #NOTE 变电站主变出线出力时间序列
    powerSubstationFileResult=open(resultPowerSubstationPath,"w")
    #NOTE 负荷供电状况时间序列
    loadEnergrizedFileResult=open(resultLoadEnergrized,"w")

    println("!!!阶段 "* timeCell*"结束************")
    println(fileResult,"!!!阶段 "* timeCell*"结束************")
    println("!!!以下是阶段 "* timeCell*"报告************")
    println(fileResult,"!!!以下是阶段 "* timeCell*"报告************")
    println("!!!阶段 "* timeCell*"开始时间=$strCurrentStartTime 结束时间=$strCurrentEndTime 共计耗时=$strRunTime  ************")
    println(fileResult,"!!!阶段 "* timeCell*"开始时间=$strCurrentStartTime 结束时间=$strCurrentEndTime 共计耗时=$strRunTime  ************")

    println("!!!阶段 "* timeCell* "优化最终状态标志位为：$state")
    println(fileResult,"!!!阶段 "* timeCell* "优化最终状态标志位为：$state")
    println("!!!阶段 "* timeCell* "最优值:",objective_value(bus33Reconfiguration))
    println(fileResult,"!!!阶段 "* timeCell* "最优值:",objective_value(bus33Reconfiguration))

    #输出初始开关状态
    println("!!!阶段 "* timeCell*"  初始各线路开关状态为: ")
    println(fileResult,"!!!阶段 "* timeCell*"  初始各线路开关状态为: ")
    n_initStateSwitch=1
    for line in lines
        println("****!!!阶段 "* timeCell*"  线路 $line 的初始状态为：",stateInit[n_initStateSwitch])
        println(fileResult,"****!!!阶段 "* timeCell*"  线路 $line 的初始状态为：",stateInit[n_initStateSwitch])
        n_initStateSwitch+=1
    end
    

    #输出各时刻各开关状态
    #横向节点\线路  纵向时间 下同 
    switchOperationPrintToFile=Array{Int64}(undef,points,length(lines))
    println("!!!以下是阶段 "*timeCell*"  各规划时刻各线路开关状态：")
    println(fileResult,"!!!以下是阶段 "*timeCell*"  各规划时刻各线路开关状态：")
    for time in 1:points
        n=0
        for line in lines
            n+=1
            if checkLineAlive(line,lineData,time,lineFLAG)
                switchOperationPrintToFile[time,n]=round(Int64,(value(alpha[(line[1],line[2],time)])))
                println("****阶段 "*timeCell*"线路 $line 在 $time 规划时刻的状态为：",switchOperationPrintToFile[time,n])
                println(fileResult,"****阶段 "*timeCell*"线路 $line 在 $time 规划时刻的状态为：",switchOperationPrintToFile[time,n])
            else
                #ATTITION! 若线路当前不可用，状态置为-999
                println("****阶段 "*timeCell*"线路 $line 在 $time 规划时刻故障，已断开！")
                println(fileResult,"****阶段 "*timeCell*"线路 $line 在 $time 规划时刻故障，已断开！")
                switchOperationPrintToFile[time,n]=-999
            end
        end        
    end
    array2dPrint2File(switchOperationPrintToFile,switchOperationFileResult)
    
    println("!!!以下是阶段 "*timeCell*" 各规划时刻发生动作的线路开关统计：")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 各规划时刻发生动作的线路开关统计：")
    for time in 1:points
        ns=0
        for line in lines
            #开关动作= lambda值变为1，且 不是由于上一步1状态被故障置为0的情况
            #if round(Int64,value(lambda[(line[1],line[2],time)]))==1
            if linesFault==()
                if round(Int64,value(lambda[(line[1],line[2],time)]))==1 
                    ns+=1
                    println("****阶段 "*timeCell*" 线路 $line 上的断路器在 $time 规划时刻进行了1次操作")
                    println(fileResult,"****阶段 "*timeCell*" 线路 $line 上的断路器在 $time 规划时刻进行了1次操作")
                end
            else
                if round(Int64,value(lambda[(line[1],line[2],time)]))==1 & !(round(Int64,value(alpha[(line[1],line[2],time-1)]))==1 & ((line[1],line[2],time) in linesFault))
                    ns+=1
                    println("****阶段 "*timeCell*" 线路 $line 上的断路器在 $time 规划时刻进行了1次操作")
                    println(fileResult,"****阶段 "*timeCell*" 线路 $line 上的断路器在 $time 规划时刻进行了1次操作")    
                end
            end
        end
        if ns==0
            println("****注意！！阶段 "*timeCell*" 配电网在 $time 规划时刻没有进行任何开关操作")
        end
    end
    

    #输出各节点电压幅值
    nodeVoltagePrintToFile=Array{Float64}(undef,points,length(nodes))
    println("!!!以下是阶段 "*timeCell*" 各规划时刻各节点电压幅值的标幺值：")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 各规划时刻各节点电压幅值的标幺值：")
    for time in 1:points
        for node in nodes
            nodeVoltagePrintToFile[time,node]=sqrt(value(nodeSquVoltage[(node,time)]))
            println("****阶段 "*timeCell*" 在$time 规划时刻，节点 $node 的电压幅值(pu.)为：",nodeVoltagePrintToFile[time,node])
            println(fileResult,"****阶段 "*timeCell*" 在$time 规划时刻，****节点 $node 的电压幅值(pu.)为：",nodeVoltagePrintToFile[time,node])
        end
    end
    array2dPrint2File(nodeVoltagePrintToFile,nodeVoltageFileResult)

    #输出各节点注入功率
    println("!!!以下是阶段 "*timeCell*" 各规划时刻各节注入有功功率的真实值：")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 各规划时刻各节注入有功功率的真实值：")
    for time in 1:points
        for node in nodes
            println("****阶段 "*timeCell*" 在$time 规划时刻，节点 $node 的注入有功功率(kW)为：",base_S/1000*value(injectionActivePower[(node,time)]))
            println(fileResult,"****阶段 "*timeCell*" 在$time 规划时刻，节点 $node 的注入有功功率(kW)为：",base_S/1000*value(injectionActivePower[(node,time)]))
        end
    end
    println("!!!以下是阶段 "*timeCell*" 各规划时刻各节注入无功功率的真实值：")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 各规划时刻各节注入无功功率的真实值：")
    for time in 1:points
        for node in nodes
            println("****阶段 "*timeCell*" 在$time 规划时刻，节点 $node 的注入无功功率(kVar)为：",base_S/1000*value(injectionReactivePower[(node,time)]))
            println(fileResult,"****阶段 "*timeCell*" 在$time 规划时刻，节点 $node 的注入无功功率(kVar)为：",base_S/1000*value(injectionReactivePower[(node,time)]))
        end
    end

    #输出各线路电流幅值
    lineCurrentPrintToFile=Array{Float64}(undef,points,length(lines))
    println("!!!以下是阶段 "*timeCell*" 各规划时刻各线路通过电流的真实值：")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 各规划时刻各线路通过电流的真实值：")
    for time in 1:points
        i=0
        for line in lines
            i+=1
            #开方
            lineCurrentPrintToFile[time,i]=base_I*sqrt(abs(value(lineSquCurrent[(line[1],line[2],time)])))    
            println("****阶段 "*timeCell*" 在$time 规划时刻，线路 $line 上流过的电流大小(A)为：",lineCurrentPrintToFile[time,i])
            println(fileResult,"****阶段 "*timeCell*" 在$time 规划时刻， 线路 $line 上流过的电流大小(A)为：",lineCurrentPrintToFile[time,i])
        end
    end
    array2dPrint2File(lineCurrentPrintToFile,lineCurrentFileResult)

    lineLossPrintToFile=Array{Float64}(undef,points,length(lines))
    println("!!!以下是阶段 "*timeCell*" 各规划时刻网损真实值统计情况：")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 各规划时刻网损真实值统计情况：")
    for time in 1:points 
        i=0
        loss_T_System=0
        for line in lines
            i+=1
            lineLossPrintToFile[time,i]=base_S/1000*rLineDict[line]*value(lineSquCurrent[(line[1],line[2],time)])
            println("****在阶段 "*timeCell* " 线路 $line 在 $time 规划时刻的网损大小(kW)为：",lineLossPrintToFile[time,i])
            println(fileResult,"****在阶段 "*timeCell* " 线路 $line 在 $time 规划时刻的网损大小(kW)为：",lineLossPrintToFile[time,i])
            loss_T_System+=base_S*rLineDict[line]/1000*value(lineSquCurrent[(line[1],line[2],time)])
        end
        println("****在阶段"*timeCell * "$time 规划时刻的总体网损大小(kW)为：",loss_T_System)    
        println(fileResult,"****在阶段"*timeCell * "$time 规划时刻的总体网损大小(kW)为：",loss_T_System)    
    end
    array2dPrint2File(lineLossPrintToFile,lineLossFileResult)

    if tsMTalive!=()
        #ATTITION!列号仅代表listMT中元素的序号
        powerMTPrintToFile=Array{Complex}(undef,points,length(listMT))
        #首先置零
        fill!(powerMTPrintToFile,0+0im)
        println("!!!以下是阶段 "*timeCell*" 各MT在各规划时刻发出功率真实值统计情况：")
        println(fileResult,"!!!以下是阶段 "*timeCell*" 各MT在各规划时刻发出功率真实值统计情况：")
        for it in tsMTalive
            powerMTPrintToFile[it[2],mtNode2AliveDataColumn(it[1],DataPath)]=complex(base_S/1000*(value(activePowerMT[(it[1],it[2])]))
                                        ,base_S/1000*(value(reactivePowerMT[(it[1],it[2])])))
            println("****阶段 "*timeCell*" MT $(it[1]) 在 $(it[2]) 规划时刻的发出有功大小(kW)为：",base_S/1000*(value(activePowerMT[(it[1],it[2])])))
            println("****阶段 "*timeCell*" MT $(it[1]) 在 $(it[2]) 规划时刻的运行成本为：",mtInfDict[it[1]][5]*(value(activePowerMT[(it[1],it[2])])))
            println(fileResult,"****阶段 "*timeCell*" MT $(it[1]) 在 $(it[2]) 规划时刻的发出有功大小(kW)为：",base_S/1000*(value(activePowerMT[(it[1],it[2])])))
            println(fileResult,"****阶段 "*timeCell*" MT $(it[1]) 在 $(it[2]) 规划时刻的运行成本为：",mtInfDict[it[1]][5]*(value(activePowerMT[(it[1],it[2])])))
            println("****阶段 "*timeCell*" MT $(it[1]) 在 $(it[2]) 规划时刻的发出无功大小(kVar)为：",base_S/1000*(value(reactivePowerMT[(it[1],it[2])])))
            println(fileResult,"****阶段 "*timeCell*" MT $(it[1]) 在 $(it[2]) 规划时刻的发出无功大小(kVar)为：",base_S/1000*(value(reactivePowerMT[(it[1],it[2])])))
        end
        array2dPrint2File(powerMTPrintToFile,powerMTFileResult)
    end

    powerSubstationPrintToFile=Array{Complex}(undef,points,length(listSub))
    #ATTITION!列号仅代表listSub中元素的序号
    println("!!!以下是阶段 "*timeCell*" 变电站主变出线功率情况")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 变电站主变出线功率情况")
    for time in 1:points
        i=0
        for subNode in listSub
            i+=1
            powerSubstationPrintToFile[time,i]=Complex(base_S/1000*(value(injectionActivePower[(subNode,time)])),base_S/1000*(value(injectionReactivePower[(subNode,time)])))
            println("****阶段 "*timeCell*" 变电站 $subNode 在 $time 规划时刻从上级电网输入的有功大小(kW)为：",real(powerSubstationPrintToFile[time,i]))
            println(fileResult,"****阶段 "*timeCell*" 变电站 $subNode 在 $time 规划时刻从上级电网输入的有功大小(kW)为：",real(powerSubstationPrintToFile[time,i]))
            println("****阶段 "*timeCell*" 变电站 $subNode 在 $time 规划时刻从上级电网输入的无功大小(kVar)为：",imag(powerSubstationPrintToFile[time,i]))
            println(fileResult,"****阶段 "*timeCell*" 变电站 $subNode 在 $time 规划时刻从上级电网输入的无功大小(kVar)为：",imag(powerSubstationPrintToFile[time,i]))
        end
    end
    array2dPrint2File(powerSubstationPrintToFile,powerSubstationFileResult)

    loadEnergrizedPrintToFile=Array{Float64}(undef,points,numNode)
    println("!!!以下是阶段 "*timeCell*" 负荷供电状态")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 负荷供电情况")
    for time in 1:points
        #根节点为-999
        for node in listSub
            loadEnergrizedPrintToFile[time,node]=-999
        end
        for node in rootFreeNodes
            loadEnergrizedPrintToFile[time,node]=value(loadRate[(node,time)])
            if loadEnergrizedPrintToFile[time,node]<0.5
                println("****警告！ 在阶段 "*timeCell*" 中 "*" $time 时刻 配电网 $node 号节点负荷削减超过50%！")
                println(fileResult, "****警告！ 在阶段 "*timeCell*" 中 "*" $time 时刻 配电网 $node 号节点负荷削减超过50%！")
            end
            println("****配电网在阶段 "*timeCell*" $time 时刻 $node 号节点负荷率为： ",loadEnergrizedPrintToFile[time,node])
            println(fileResult, "****配电网在阶段 "*timeCell*" $time 时刻 $node 号节点负荷率为： ",loadEnergrizedPrintToFile[time,node])
        end    
    end
    array2dPrint2File(loadEnergrizedPrintToFile,loadEnergrizedFileResult)

    println("!!!以下是阶段 "*timeCell*" 负荷停电损失情况")
    println(fileResult,"!!!以下是阶段 "*timeCell*" 负荷停电损失情况")
    for t in 1:points
        sumBlackoutLoss=0
        for node in nodes
            blackoutLoss=failureLoss*consumeP[t,node]*(1-value(loadRate[(node,t)]))
            sumBlackoutLoss+=blackoutLoss
            println("****配电网在阶段 "*timeCell*" $t 规划时刻 节点 $node 负荷停电损失为： ",blackoutLoss)
            println(fileResult,"****配电网在阶段 "*timeCell*" $t 规划时刻 节点 $node 负荷停电损失为： ", blackoutLoss)
        end
        println(fileResult,"****配电网在阶段 "*timeCell*" 总负荷停电损失为： ", sumBlackoutLoss)
    end
    
    println("###########阶段 " *timeCell* "#######结束##########################")
    println(fileResult,bus33Reconfiguration)
    println("###########阶段 " *timeCell* " 报告保存完成##########################")

    # println("******#测试用！！！ 显示F的情况***")
    # println(fileResult,"******#测试用！！！ 显示F的情况***")
    # for ijt in linesAlive
    #     println("****在 $(ijt[3]) 规划时刻 线路$((ijt[1],ijt[2])) 的F值为： ",round(Int64,value(apparentFictitiousFlow[(ijt[1],ijt[2],ijt[3])])))
    #     println(fileResult,"****在 $(ijt[3]) 规划时刻 线路$((ijt[1],ijt[2])) 的F值为： ",round(Int64,value(apparentFictitiousFlow[(ijt[1],ijt[2],ijt[3])])))
    # # end  
    # for t in 0:points
    #     for ij in lines
    #         println("线路 $(ij[1]) ,$(ij[2]) 在 $t 时的alpha值为 ",value.(alpha[(ij[1],ij[2],t)]))
    #     end
    # end

    # for t in 1:points
    #     for ij in lines
    #         println("线路 $(ij[1]) ,$(ij[2]) 在 $t 时的lambda值为 ",value.(lambda[(ij[1],ij[2],t)]))
    #     end
    # end
   
   
    #关闭文件读写    
    close(fileResult)
    close(switchOperationFileResult)
    close(powerMTFileResult)
    close(nodeVoltageFileResult)
    close(lineCurrentFileResult)
    close(lineLossFileResult)
    close(powerSubstationFileResult)

    #将最后一步规划的开关状态更新stateInit
    updateStateInit(switchOperationPrintToFile,points,stateInit)
end