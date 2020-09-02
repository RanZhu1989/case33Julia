#NOTE  2020.8.28 v0.5 
#NOTE  Author：Ran Zhu @ School of Cyber Engineering and Science, SEU
#NOTE  此.jl脚本用于测试评估锥规划模型下正常运行和一般故障下的自适应三相平衡配网重构单步规划和长时间尺度动态规划数据生成器的性能
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


import Pkg
import MathOptInterface
using DelimitedFiles
using JuMP,MosekTools,Convex

####################静态数据集路径###########################


#ATTITION!现在填入的是测试用单步规划数据集
#ATTITION!注意，原始数据需要在.jl脚本中转换为标幺值，后面计算采用标幺值计算

const loadFilePath="33BUSdata//forTEST//TEST_ts_33bus_load.csv"
const pvFilePath="33BUSdata//forTEST//TEST_ts_33bus_PV_MOD.csv"
const paraLinePath="33BUSdata//forTEST//TEST_para_33bus_line_Li2019TSG.csv"
const paraNodePath="33BUSdata//forTEST//TEST_para_33bus_node_Li2019TSG.csv"


#####################报告文件路径##############################
#ATTITION! 必须手动创建report目录  Julia下的文件操作尚不完善
#ATTITION! mkdir若目录已存在则会报错
#mkdir(".//report//")
const reportPath=".//report//"
const resultPath="result重构.txt"
const resultSwitchOpearingPath=reportPath * "开关动作_" * resultPath
const resultPowerMTPath=reportPath *"MT出力_kW_" * resultPath
const resultNodeVoltagePath=reportPath * "节点电压_pu._" * resultPath
const resultLineCurrentPath=reportPath * "线路电流_A_" * resultPath
const resultLossPath=reportPath * "线路网损_" * resultPath
const resultPowerSubstationPath=reportPath * "变电站主变出线出力_" *resultPath
const resultLoadEnergrized=reportPath * "负荷供电状态_" *resultPath

#######################全局常量##############################
#基准容量
#用于计算标幺值
#Z_B=(U_B)^2 / S_B
#S_B = sqrt(3) * U_B * I_B
#r* =r/Z_B
#x* =x/Z_B
#p* =p/S_B
#q* =q/S_B

#NOTE: 设置基准容量 => 10kV等级配网设置为100MVA
const base_S=100e6
#NOTE: 设置基准电压 => 依据具体算例设置
const base_V=12.66e3
const base_I=(base_S)/sqrt(3)/base_V
const base_Z=(base_I)^2/base_S

#功率因数
#TODO：计算真实的大数据序列时有用
const powerFactor=0.95

#规划时间尺度 26周
#points=4368
#单步规划 points=1
const points=1
#节点负荷时间序列矩阵 P or Q(i,t) 行i=节点号 列t=时间
#FIXME补全读取函数
#ATTITION!测试时不要运行下面这条语句
#const consumeP,consumeQ=readConsumePower(loadFilePath,powerFactor,points)
#PV有功出力时间序列矩阵 P(i,t) 行i=节点号 列t=时间
#ATTITION!测试时不要运行下面这条语句
#FIXME补全读取函数
#const PvP=readPVPower(pvFilePath,points)

#网损电价 单位monemy / kW*h 千瓦时
#NOTE: 来源 2020, Gao ,TSG
const priceLoss=0.13/(1000/base_S)#FIXME:网损成本需要设置，找参考文献！
#NOTE: 来源：？？？
#购电电价 单位monemy / kW*h 千瓦时
const priceGrid=0.13/(1000/base_S)#FIXME:购电电价需要设置，找参考文献！
#操作损耗
const costSwitch=0.1#FIXME:操作成本需要设置，找参考文献！
#停电单位损失
const failureLoss=10/(1000/base_S)#FIXME:停电成本需要设置，找参考文献！
const unitLossPenaltyCoefficient=1000#FIXME:停电成本需要设置，找参考文献！
#节点电压上下界,原文中电压等级为12kV
#ATTITION!为线性化采用了电压幅值的平方作为变量！
const lowSquVNode=(0.90)^2
const highSquVNode=(1.1)^2
#变电站主变出线电压幅值
const voltageSquSub=(1)^2
#具备黑启动能力的DG的电压幅值
const voltageBlackStartDG=1.00
#ATTITION!每次开关最大动作次数 动态规划时使用这一条语句
#const maxNS=2
#ATTITION!无限制单步的优化用下面这条语句
#NOTE 设置为0 csv文件中修改上下限约束 可以作为潮流计算用
const maxNS=99
#大M 至少大于节点数即可
const bigM=50
####################ATTITION!以下为形成测试数据的语句，生成.jl文件时注释掉
#读取TESTload
#ATTITION!.csv文件填入真实值，读取时转换为标幺值
#已测试
readTESTLoadData=readdlm(loadFilePath,',',header=true)
TESTLoadData=readTESTLoadData[1]
TESTlistP=TESTLoadData[:,3]
TESTlistQ=TESTLoadData[:,4]
global consumeP=zeros(33,points)
global consumeQ=zeros(33,points)
for iTESTLoad in 1:33
    for tTESTLoad in 1:points
        consumeP[iTESTLoad,tTESTLoad]=(TESTlistP[iTESTLoad])/base_S
        consumeQ[iTESTLoad,tTESTLoad]=(TESTlistQ[iTESTLoad])/base_S    
    end
end
#读取TESTPV
#已测试
readTESTPVData=readdlm(pvFilePath,',',header=true)
TESTPVData=readTESTPVData[1]
global PVP=zeros(33,points)
TESTlistPVP=TESTPVData[:,3]
for iTESTPV in 1:33
    for tTESTLoad in 1:points
        PVP[iTESTPV,tTESTLoad]=(TESTlistPVP[iTESTPV])/base_S
    end
end

##############################功能函数#############################################
#= function readConsumePower(filePath,powerFactor,points)
    #TODO：补全功能
    # 在输入时间尺度以及功率因数下返回各节点有功负载、无功负载.
    #格式：
    #
    return P,Q
end =#

#= function readPVPower(filePath,points)
    #TODO：补全功能
    # 在输入时间尺度以及功率因数下返回PV有功输出(假设该PV采用功率因数控制，不发无功).
    return P
end =#

#读取初始拓扑配置
#返回存活的常闭开关元组，常开开关元组
function readtopologyData(data)
    listIC=[]
    listIO=[]
    #常闭
    indexIC=findall(x->x==1,data[:,9])
    #常开
    indexIO=findall(x->x==0,data[:,9])
    for i in 1:length(indexIC)
        #检测存活性
        if round(Int64,data[indexIC[i],8])==1
            push!(listIC,(round(Int64,data[indexIC[i],2]),round(Int64,data[indexIC[i],3])))
        end
    end
    for i in 1:length(indexIO)
        if round(Int64,data[indexIO[i],8])==1
            push!(listIO,(round(Int64,data[indexIO[i],2]),round(Int64,data[indexIO[i],3])))
        end
    end
    return Tuple(listIC), Tuple(listIO)
end

#将二维数组打印到文件
function array2dPrintToFile(printArray,fileOpened)
    for i in 1:(size(printArray))[1]
        for j in 1:(size(printArray))[2]
            print(fileOpened,printArray[i,j])
            print(fileOpened," ")
        end
        println(fileOpened)
    end
end

function linesijTOlinesIJJI(lines,lineData)
    #创建一个将(ij)并(ji)的元组集合
    result=[]
    for line in lines
        if checkLineAlive(line,lineData)
            #取元组中数据
            temp1=(line[1],line[2])
            temp2=(line[2],line[1])
            push!(result,temp1)
            push!(result,temp2)
        end    
    end
    return Tuple(result)
end

function checkNodeFlagAlive(i,Data,j)
    #已测试
    #检查一个节点是否存活 返回布尔量
    #存活指负荷\发电机是否退出
    #Data是node.csv文件，j=10 负载，j=11 PV j=12 MT
    if Data[i,j]==1
        return true
    else
        return false
    end
end

function checkLineAlive(line,Data)
    #已测试
    #检查一条线路(i,j)是否存活
    for i in 1:length(Data[:,1])
        if (Data[i,2],Data[i,3])==line
            if Data[i,8]==1
                return true
            end        
        end
    end
    return false
end

function findSub(Data)
    #已测试
    #寻找并返回变电站节点编号的元组
    find=[]
    i=0
    for j in Data[:,2]
        i+=1
        if round.(Int64,j)==1
           push!(find,i) 
        end
    end
    return Tuple(find)
end

function findPV(Data)
    #已测试
    #寻找并返回PV节点编号的元组
    find=[]
    i=0
    for j in Data[:,3]
        i+=1
        if round.(Int64,j)==1
           push!(find,i) 
        end
    end
    return Tuple(find)
end

function findPureLoad(Data,numNode)
    #已测试
    #寻找并返回纯负载节点编号的元组
    Sub=Set(findSub(Data))
    PV=Set(findPV(Data))
    MT=Set(keys(findMT(Data,base_S)))
    #并集运算
    temp=union(Sub,PV)
    temp=union(temp,MT)
    ALL=Set(i for i in 1:numNode)
    #差集运算
    Load=setdiff(ALL,temp)
    return Tuple(Load)
end

function findMT(Data,base_S)
    #已测试
    #寻找并返回:
    #生成MT编号、运行成本、MT出力上下界的字典
    #ATTITION!生成的字典将实际参数值转换为标幺值
    findNode=[]
    findLP=[]
    findHP=[]
    findLQ=[]
    findHQ=[]
    findCost=[]
    findAlive=[]
    i=0
    for j in Data[:,4]
        i+=1
        if round.(Int64,j)==1
           push!(findNode,i)
           push!(findLP,(Data[i,5])/base_S) 
           push!(findHP,(Data[i,6])/base_S) 
           push!(findLQ,(Data[i,7])/base_S) 
           push!(findHQ,(Data[i,8])/base_S) 
           push!(findAlive,Data[i,12]) 
           #ATTITION!在这里转换为标幺值下的价格
           push!(findCost,(Data[i,9])/(1000/base_S)) 
        end
    end
    TfindNode=Tuple(findNode)
    MTdata=Tuple((findLP[k],findHP[k],findLQ[k],findHQ[k],findCost[k],findAlive[k]) for k in 1:length(TfindNode))
    MTdict=Dict(TfindNode .=>MTdata)
    return MTdict
end

function linesToTimePair(lines,points,start)
    #已测试
    #将(i,j)组成的元组转换为(i,j,t)组成的大型元组
    #ATTITION!这里时间可以设置从0或1开始
    result=[]
    for pair in lines
        for t in start:points
            #取元组中数据
            temp=(pair[1],pair[2],t)
            push!(result,temp)
        end
    end
    ijtPairTuple=Tuple(result)
    return ijtPairTuple
end

function ijjiTimePair(lines,points,start)
    #已测试
    #将(i,j)组成的元组转换为(i,j,t),(j,i,t)组成的大型元组
    #ATTITION!这里时间可以设置从0或1开始
    result=[]
    for pair in lines
        for t in start:points
            #取元组中数据
            temp1=(pair[1],pair[2],t)
            temp2=(pair[2],pair[1],t)
            push!(result,temp1)
            push!(result,temp2)
        end
    end
    ijjitPairTuple=Tuple(result)
    return ijjitPairTuple
end

function time0ijtPair(lines)
    #已测试
    #生成(i,j,0)，
    result=[]
    for pair in lines
        temp=(pair[1],pair[2],0)
        push!(result,temp)
    end
    return Tuple(result)    
end

function ij0SET(Data)
    #已测试
    #分别返回联络开关的(i,j,0) 以及普通线路的(i,j,0)
    resultTie=[]
    resultCom=[]
    i=0
    for j in Data[:,7] 
        i+=1
        if round.(Int64,j)==1
            startNode=round.(Int64,Data[i,2])
            endNode=round.(Int64,Data[i,3])
            push!(resultTie,(startNode,endNode,0))
        else
            startNode=round.(Int64,Data[i,2])
            endNode=round.(Int64,Data[i,3])
            push!(resultCom,(startNode,endNode,0))
        end
    end
    return Tuple(resultTie),Tuple(resultCom)
end

function listPoints(points)
    #已测试
    #返回一个1~points的元组
    temp=[i for i in 1:points]
    return Tuple(temp)
end

function time1itPair(numNode,points)
    #已测试
    #返回一个(i，t)组成的元组
    result=[]
    for i in 1:numNode
        for t in 1:points
            temp=(i,t)
            push!(result,temp)
        end
    end
    return Tuple(result)
end

function findijNeighborNode(i,lines,Data)
    #已测试
    #返回节点i与所有可用邻居组成的(i,j)元组
    #这里的可用指：相连的线路存活为真
    #i=节点标号 lines=(i,j)元组集合 Data=线路
    result=[]
    for line in lines
        if checkLineAlive(line,Data)
            if i in line
                push!(result,line)
            end
        end
    end
    return Tuple(result)
end

function rootFreeFindijtNeighborNode(i,lines,points)
    #已测试
    #ATTITION!这是针对(i,j)U(j,i)的一个“向后"的搜索
    #NOTE 某些情况下可能存在没有邻居的情况
    #返回节点i与所有可用邻居组成的(i,j,t)元组
    #这里的可用指：相连的线路存活为真
    #i=节点标号 lines=(i,j)元组集合 Data=读取的线路数据
    result=[]
    for line in lines
        if i==line[1]
            for t in 1:points
                temp=(line[1],line[2],t)
                push!(result,temp)
            end
        end
    end
    if result!=[]
        return Tuple(result)
    else
        return ()
    end
end

function findjkForwardNeighborPair(j,lines,Data)
    #已测试
    #返回节点j在给定方向与下游邻居节点组成的(j,k)元组
    #方向实际由csv文件的起讫点给出
    result=[]
    allPair=findijNeighborNode(j,lines,Data)
    for pair in allPair
        if pair[1]==j
            push!(result,pair)
        end
    end
    return Tuple(result)
end

function findijBackforwardNeighborPair(j,lines,Data)
    #已测试
    #返回节点j在给定方向与上游邻居节点组成的(i,j)元组
    #方向实际由csv文件的起讫点给出
    result=[]
    allPair=findijNeighborNode(j,lines,Data)
    for pair in allPair
        if pair[2]==j
            push!(result,pair)
        end
    end
    return Tuple(result)
end

function findTupleAliveLines(lines,lineData)
    #返回存活的线路组成的(i,j)元组
    result=[]
    for line in lines
        if  checkLineAlive(line,lineData)
            push!(result,line)
        end
    end
    return  Tuple(result)
end

function findTupleFaultLines(lines,lineData)
    #返回故障的线路组成的(i,j)元组
    result=[]
    for line in lines
        if  checkLineAlive(line,lineData)==false
            push!(result,line)
        end
    end
    return  Tuple(result)
end

function nodeMTAliveitPair(listMT,dic)
    result=[]
    for node in listMT
        if dic[node][6]==1
            for t in 1:points
                push!(result,(node,t))
            end
        end
    end
    return  Tuple(result)
end
###############################

function dpSolverReconfiguraiton33Bus()
    
    ####################读取必要数据文件########################################

    ##TODO 下面两个路径可作为参数传入
    #文件路径
    #paraLinePath="33BUSdata//para_33bus_line.csv"
    #paraNodePath="33BUSdata//para_33bus_node.csv" 

    #ATTITION!仅供测试用！形成脚本时注意注释掉
    #读取线路参数文件
    #reference:  M.E.Baran,1989,TPS
    readLineData=readdlm(paraLinePath,',',header=true)
    lineData=readLineData[1]
    listLine=round.(Int64,lineData[:,1])
    startNode=round.(Int64,lineData[:,2])
    endNode=round.(Int64,lineData[:,3])
    stateInit=round.(Int64,lineData[:,9])
    lineResistance=(lineData[:,4])/base_Z
    lineImpedance=(lineData[:,5])/base_Z
    lineSquCurrentMAX=(lineData[:,6])/base_I
    #读取联络开关位置，备用
    listFlagTieSwitch=lineData[:,7]
    #生成线路相关数据结构
    #读取待求解问题常闭开关，常开开关
    listIC,listIO=readtopologyData(lineData)
    #总节点个数，线路条数
    numNode=max(maximum(startNode),maximum(endNode))
    numLine=length(listLine)
    #遍历系统中所有节点的迭代器
    nodes=1:numNode
    #遍历非根节点的迭代器
    #ATTITION!：这里人工指定了root=1
    rootFreeNodes=2:numNode
    #(i,t) 节点-时间 二元组
    itPair=time1itPair(numNode,points)
    #生成标准情况下的线路(i,j)元组
    lines=Tuple((startNode[i],endNode[i]) for i in 1:numLine)
    #生成存活的线路(i,j)元组
    linesAlive=findTupleAliveLines(lines,lineData)
    #生成故障线路(i,j)元组
    linesFault=findTupleFaultLines(lines,lineData)
    #标准情况下逆序的线路(j,i)元组
    reverseLines=Tuple((endNode[i],startNode[i]) for i in 1:numLine)
    #考虑线路生存状态的(i,j)U(j,i)
    ijjiLines=linesijTOlinesIJJI(lines,lineData)
    #线路电阻值与阻抗的字典
    rLineDict=Dict(lines .=>lineResistance)
    xLineDict=Dict(lines .=>lineImpedance)
    #线路最大载流量的字典
    maxCurrentLineDict=Dict(lines .=>lineSquCurrentMAX)
    #常闭开关所在线路的 (i,j,t)三元组穷举 从t=1开始
    time1ICSwitchPair=linesToTimePair(listIC,points,1)
    #常开开关所在线路的(i,j,t)三元组穷举 从t=1开始
    time1ICOwitchPair=linesToTimePair(listIO,points,1)
    #存活的(i,j,t)三元组穷举，i,j∈line，t=0或1:points
    time0NodePairAlive=linesToTimePair(linesAlive,points,0)
    time1NodePairAlive=linesToTimePair(linesAlive,points,1)
    #所有的(i,j,t)三元组穷举
    time1NodePair=linesToTimePair(lines,points,1)
    #(i,j,t) (j,i,t)三元组 t=1：points
    ijtjit1Pair=ijjiTimePair(lines,points,1)
    #(i,j,0)三元组
    ij0Pair=time0ijtPair(linesAlive)
    #联络开关(i,j,0) 以及普通支路(i,j,0) 没用
    ij0TiePair,ij0ComPair=ij0SET(lineData)
    #读取节点参数文件
    readNodeData=readdlm(paraNodePath,',',header=true)
    nodeData=readNodeData[1]
    #找出变电站节点元组
    listSub=findSub(nodeData)
    #找出PV节点元组
    listPV=findPV(nodeData)
    #找出并创建MT节点信息的字典
    #字典说明：键值=节点编号；值=(有功下界，有功上界，无功下界，无功上界，成本C, 生存)
    mtInfDict=findMT(nodeData,base_S)
    #找出MT节点元组
    listMT=Tuple(Set(keys(mtInfDict)))
    #有效MT节点的(i,t)
    itPairMT=nodeMTAliveitPair(listMT,mtInfDict)
    #找出普通负载节点元组
    listPureLoad=findPureLoad(nodeData,numNode)

    ###################建立JuMP模型##############################
    bus33Reconfiguration=Model(with_optimizer(Mosek.Optimizer))
  
    ######################变量设置##################################
   
    #表示线路通电状态的lineEnergized 二进制矩阵 实际意思是线路闭合且通电 1=线路闭合带电 0=线路闭合但不带电通电
    @variable(bus33Reconfiguration,lineEnergized[ijt in time1NodePair],Bin,base_name="闭合线路是否带电lineEnergized")

    #表示节点通电状态的nodeEnergized 二进制矩阵
    @variable(bus33Reconfiguration,nodeEnergized[it in itPair],Bin,base_name="节点是否带电nodeEnergized")

    #表示线路断路器即使闭合还是不带电的closedLineDeEnergized 二进制矩阵 1=线路闭合但不带电 0=线路闭合且通电
    #ATTITION! 只针对初始状态为闭合的开关有效
    #NOTE 在目标最大化负载生存率，以及最小化操作数的条件下，大部分情况下以下逻辑成立：
    #NOTE 闭合开关=>带电  打开开关=>失电
    #NOTE 但对于初态是（闭合，带电）的状态，有可能变为（闭合，不带电） 见 2019, Li, TSG
    @variable(bus33Reconfiguration,closedLineDeEnergized[ijt in time1NodePair],Bin,base_name="闭合线路是否失电closedLineDeEnergized")

    #用于构建辐射状拓扑约束的中ST约束辅助变量b_{ijt}
    #ATTITION!bAuxiliary包含(i,j)和对应元素交换后的(j,i)
    @variable(bus33Reconfiguration,bAuxiliary[ijt in ijtjit1Pair],Bin,base_name="辅助变量b")

    #MT有功出力P_mt_{i t}
    @variable(bus33Reconfiguration,activePowerMT[it in itPairMT],base_name="MT有功出力")
    #MT无功出力Q_mt_{i t}
    @variable(bus33Reconfiguration,reactivePowerMT[it in itPairMT],base_name="MT有功出力")

    #注入有功\无功功率 P Q _{i t}
    @variable(bus33Reconfiguration,injectionActivePower[it in itPair],base_name="注入有功")
    @variable(bus33Reconfiguration,injectionReactivePower[it in itPair],base_name="注入无功")

    #首端有功\无功功率 P Q _{i j t}
    @variable(bus33Reconfiguration,apparentActivePower[ijt in time1NodePairAlive],base_name="首端有功")
    @variable(bus33Reconfiguration,apparentReactivePower[ijt in time1NodePairAlive],base_name="首端无功")

    #ATTITION! l_{i j t} 支路电流平方
    #每条支路最大载流量需要查lineData
    @variable(bus33Reconfiguration,lineSquCurrent[ijt in time1NodePairAlive],base_name="线路电流幅值平方")

    #ATTITION!节点电压平方幅值 V_{i t}
    @variable(bus33Reconfiguration,nodeSquVoltage[it in itPair],base_name="节点电压幅值平方")

    #######################补充的变量取值范围#########################

    #存在故障线路情况下的lineEnergized与closedLineDeEnergized初值设置
    for t in points
        for line in linesFault
            @constraint(bus33Reconfiguration,lineEnergized[(line[1],line[2],t)]==0)  
            @constraint(bus33Reconfiguration,closedLineDeEnergized[(line[1],line[2],t)]==0)
            @constraint(bus33Reconfiguration,bAuxiliary[(line[1],line[2],t)]==0)
            @constraint(bus33Reconfiguration,bAuxiliary[(line[2],line[1],t)]==0)
        end
    end


    #closedLineDeEnergized的定义
    #针对闭合的线路不带电的情况
    for ijt in time1ICSwitchPair
        @constraint(bus33Reconfiguration,closedLineDeEnergized[ijt]==1-nodeEnergized[(ijt[1],ijt[3])]-nodeEnergized[(ijt[2],ijt[3])]+lineEnergized[ijt]) 
    end

    #每条支路最大载流量需要查lineData
    #NOTE 加入了对于非通电线路的处理
    for line in linesAlive
        for t in 1:points
            @constraint(bus33Reconfiguration,lineSquCurrent[(line[1],line[2],t)]<=lineEnergized[(line[1],line[2],t)]*(maxCurrentLineDict[line])^2)
            @constraint(bus33Reconfiguration,lineSquCurrent[(line[1],line[2],t)]>=0)
        end
    end

    #MT出力约束
    for i in listMT
        for t in 1:points
            mtInf=mtInfDict[i]
            @constraint(bus33Reconfiguration,mtInf[1]*nodeEnergized[(i,t)]<=activePowerMT[(i,t)])
            @constraint(bus33Reconfiguration,activePowerMT[(i,t)]<=mtInf[2]*nodeEnergized[(i,t)])
            @constraint(bus33Reconfiguration,mtInf[3]*nodeEnergized[(i,t)]<=reactivePowerMT[(i,t)])
            @constraint(bus33Reconfiguration,reactivePowerMT[(i,t)]<=mtInf[4]*nodeEnergized[(i,t)])
        end
    end

    #节点电压约束
    #NOTE 加入了对于非通电线路的处理
    for it in itPair
        @constraint(bus33Reconfiguration,nodeSquVoltage[it]<=nodeEnergized[it]*highSquVNode)
        @constraint(bus33Reconfiguration,nodeSquVoltage[it]>=nodeEnergized[it]*lowSquVNode)
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
    for node in listMT
        for t in 1:points
            @constraint(bus33Reconfiguration,nodeSquVoltage[(node,t)]==voltageBlackStartDG)
        end
    end

    
    #注入功率：P_i=P发-P用
    #注意这里需要分以下情况：
    #1.纯负载节点或纯变电站节点;2.含PV的节点;3.含MT的节点
    #对于一般节点 
    #编译通过
    #NOTE 加入了对于非通电线路的处理
    for i in listPureLoad
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==-consumeP[i,t]*nodeEnergized[(i,t)])
            @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==-consumeQ[i,t]*nodeEnergized[(i,t)])
        end
    end    

    #对于含PV的节点
    if listPV!=()
        for i in listPV
            for t in 1:points
                @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==PVP[i,t]-consumeP[i,t]*nodeEnergized[(i,t)])
                @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==-consumeQ[i,t]*nodeEnergized[(i,t)])
            end
        end
    end
    
    #对于含MT的节点
    if listMT!=()
        for i in listMT
            for t in 1:points
                mtInf=mtInfDict[i]
                @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==activePowerMT[(i,t)]-consumeP[i,t]*nodeEnergized[(i,t)])
                @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==reactivePowerMT[(i,t)]-consumeQ[i,t]*nodeEnergized[(i,t)])
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
    #自带线路存活检测
    for node in nodes
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(node,t)]
            ==sum(apparentActivePower[(j,k,t)] for (j,k) in findjkForwardNeighborPair(node,lines,lineData) if j==node)
            -sum(apparentActivePower[(i,j,t)]-rLineDict[(i,j)]*lineSquCurrent[(i,j,t)] 
            for (i,j) in findijBackforwardNeighborPair(node,lines,lineData) if j==node))
        
            @constraint(bus33Reconfiguration,injectionReactivePower[(node,t)]
            ==sum(apparentReactivePower[(j,k,t)] for (j,k) in findjkForwardNeighborPair(node,lines,lineData) if j==node)
            -sum(apparentReactivePower[(i,j,t)]-xLineDict[(i,j)]*lineSquCurrent[(i,j,t)] 
            for (i,j) in findijBackforwardNeighborPair(node,lines,lineData) if j==node))
        end
    end

    #节点电压联系方程
    #编译通过
    for line in linesAlive
        for t in 1:points
            @constraint(bus33Reconfiguration,nodeSquVoltage[(line[1],t)]-nodeSquVoltage[(line[2],t)]>=-bigM*(1-lineEnergized[(line[1],line[2],t)])
            +2*(rLineDict[line]*apparentActivePower[(line[1],line[2],t)]+xLineDict[line]*apparentReactivePower[(line[1],line[2],t)])
            -((rLineDict[line])^2+(xLineDict[line]^2))*lineSquCurrent[(line[1],line[2],t)])

            @constraint(bus33Reconfiguration,nodeSquVoltage[(line[1],t)]-nodeSquVoltage[(line[2],t)]<=bigM*(1-lineEnergized[(line[1],line[2],t)])
            +2*(rLineDict[line]*apparentActivePower[(line[1],line[2],t)]+xLineDict[line]*apparentReactivePower[(line[1],line[2],t)])
            -((rLineDict[line])^2+(xLineDict[line]^2))*lineSquCurrent[(line[1],line[2],t)])
        end
    end

    #首端功率P^2+Q^2=VI^2的旋转二阶锥约束
    #编译通过
    for line in lines
        if checkLineAlive(line,lineData)
            for t in 1:points
                @constraint(bus33Reconfiguration,[lineSquCurrent[(line[1],line[2],t)],0.5*nodeSquVoltage[(line[1],t)],
                apparentActivePower[(line[1],line[2],t)],apparentReactivePower[(line[1],line[2],t)]] in RotatedSecondOrderCone())
            end
        end
    end

    #基于ST生成树约束
    #Reference: ST = J.A.Taylor,2012,TPS;
    #csv中的(start,end)中的start已被定义为end的父节点

    #ST约束
    #编译通过

    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,(bAuxiliary[(line[1],line[2],t)]+bAuxiliary[(line[2],line[1],t)])==lineEnergized[(line[1],line[2],t)])
        end
    end
    for node in rootFreeNodes
        for t in 1:points   
            #在寻找邻居节点时已自带线路存活检测功能
            @constraint(bus33Reconfiguration,sum(bAuxiliary[ijt] for ijt in rootFreeFindijtNeighborNode(node,ijjiLines,points))==nodeEnergized[(node,t)])
        end
    end

    #ATTITION!：这里假设变电站就是唯一的root节点了
    subNeiijtPair=rootFreeFindijtNeighborNode(1,ijjiLines,points)
    for ijt in subNeiijtPair
        @constraint(bus33Reconfiguration,bAuxiliary[ijt]==0)
    end
  

    #开关操作次数统计
    # Reference: 2019, Li ,TSG
    for time in points
        @constraint(bus33Reconfiguration,sum(lineEnergized[(ij[1],ij[2],time)] for ij in listIO)
        +sum(1-(lineEnergized[(ij[1],ij[2],time)]) for ij in listIC)-sum(closedLineDeEnergized[(ij[1],ij[2],time)] for ij in listIC)<=maxNS)
    end

    ###################目标函数######################################

    #编译通过
    @objective(bus33Reconfiguration,Min,sum(
    sum((mtInfDict[mtNode][5])*(injectionActivePower[(mtNode,t)]+consumeP[mtNode,t]) for mtNode in listMT)
    +sum(priceLoss*lineSquCurrent[(line[1],line[2],t)]*rLineDict[line] for line in linesAlive)
    +sum(priceGrid*injectionActivePower[(subNode,t)] for subNode in listSub)
    +sum(costSwitch*(sum(lineEnergized[(ij[1],ij[2],t)] for ij in listIO)+sum((1-lineEnergized[(ij[1],ij[2],t)]) for ij in listIC)-sum(closedLineDeEnergized[(ij[1],ij[2],t)] for ij in listIC)))
    +sum(unitLossPenaltyCoefficient*(1-nodeEnergized[(i,t)]) for i in nodes) 
    for t in 1:points))
    
    
    #运行
    optimize!(bus33Reconfiguration)

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

    println("!!!结束，优化最终状态标志位为：$state")
    println(fileResult,"!!!结束，优化最终状态标志位为：$state")
    println("!!!最优值:",objective_value(bus33Reconfiguration))
    println(fileResult,"!!!最优值:",objective_value(bus33Reconfiguration))

    #输出初始开关状态
    println("!!!初始各线路开关状态为: ")
    println(fileResult,"!!!初始各线路开关状态为: ")
    n_initStateSwitch=1
    for line in lines
        println("****线路 $line 的初始状态为：",stateInit[n_initStateSwitch])
        println(fileResult,"****线路 $line 的初始状态为：",stateInit[n_initStateSwitch])
        n_initStateSwitch+=1
    end
    

    #输出各时刻各开关状态
    switchOperationPrintToFile=Array{Int64}(undef,length(lines),points)
    println("!!!以下是各规划时刻各线路开关状态：")
    println(fileResult,"!!!以下是各规划时刻各线路开关状态：")
    for time in 1:points
        n=0
        for line in lines
            n+=1
            if checkLineAlive(line,lineData)
                if line in  listIC
                    switchOperationPrintToFile[n,time]=round(Int64,(value(lineEnergized[(line[1],line[2],time)])+value(closedLineDeEnergized[(line[1],line[2],time)])))
                    println("****线路 $line 在 $time 规划时刻的状态为：",switchOperationPrintToFile[n,time])
                    println(fileResult,"****线路 $line 在 $time 规划时刻的状态为：",switchOperationPrintToFile[n,time])
                else
                    switchOperationPrintToFile[n,time]=round(Int64,(value(lineEnergized[(line[1],line[2],time)])))
                    println("****线路 $line 在 $time 规划时刻的状态为：",switchOperationPrintToFile[n,time])
                    println(fileResult,"****线路 $line 在 $time 规划时刻的状态为：",switchOperationPrintToFile[n,time])
                end
            else
                #ATTITION! 若线路当前不可用，状态置为-999
                switchOperationPrintToFile[n,time]=-999
            end
        end    
    end
    array2dPrintToFile(switchOperationPrintToFile,switchOperationFileResult)
    
    println("!!!以下是该问题中的故障点：")
    println(fileResult,"!!!以下是该问题中故障点：")
    for line in linesFault
        println("线路 $line 发生故障")
        println(fileResult,"线路 $line 发生故障")
    end

    println("!!!以下是各规划时刻发生动作的线路开关统计：")
    println(fileResult,"!!!以下是各规划时刻发生动作的线路开关统计：")
    for time in 1:points
        ns=0
        for line in linesAlive
            if time==1
                if switchOperationPrintToFile[findall(x->x==line,lines)]!=stateInit[findall(x->x==line,lines)]
                    println("****配电网在 $time 规划时刻于线路 $line 上进行了1次开关操作")
                    ns+=1
                end
            else
                if switchOperationPrintToFile[findall(x->x==line,lines),time]!=switchOperationPrintToFile[findall(x->x==line,lines),time-1]
                    println("****配电网在 $time 规划时刻于线路 $line 上进行了1次开关操作")
                    ns+=1
                end
            end
        end
        if ns==0 
            println("****警告！配电网在 $time 规划时刻没有进行任何开关操作")
            println(fileResult,"****警告！配电网在 $time 规划时刻没有进行任何开关操作")
        end
    end
    

    #输出各节点电压幅值
    nodeVoltagePrintToFile=Array{Float64}(undef,length(nodes),points)
    println("!!!以下是各规划时刻各节点电压幅值的标幺值：")
    println(fileResult,"!!!以下是各规划时刻各节点电压幅值的标幺值：")
    for time in 1:points
        for node in nodes
            nodeVoltagePrintToFile[node,time]=sqrt(value(nodeSquVoltage[(node,time)]))
            println("****节点 $node 的电压幅值(pu.)为：",nodeVoltagePrintToFile[node,time])
            println(fileResult,"****节点 $node 的电压幅值(pu.)为：",nodeVoltagePrintToFile[node,time])
        end
    end
    array2dPrintToFile(nodeVoltagePrintToFile,nodeVoltageFileResult)

    #输出各节点注入功率
    println("!!!以下是各规划时刻各节注入有功功率的真实值：")
    println(fileResult,"!!!以下是各规划时刻各节注入有功功率的真实值：")
    for time in 1:points
        for node in nodes
            println("****节点 $node 的注入有功功率(kW)为：",base_S/1000*value(injectionActivePower[(node,time)]))
            println(fileResult,"****节点 $node 的注入有功功率(kW)为：",base_S/1000*value(injectionActivePower[(node,time)]))
        end
    end
    println("!!!以下是各规划时刻各节注入无功功率的真实值：")
    println(fileResult,"!!!以下是各规划时刻各节注入无功功率的真实值：")
    for time in 1:points
        for node in nodes
            println("****节点 $node 的注入无功功率(kVar)为：",base_S/1000*value(injectionReactivePower[(node,time)]))
            println(fileResult,"****节点 $node 的注入无功功率(kVar)为：",base_S/1000*value(injectionReactivePower[(node,time)]))
        end
    end

    #输出各线路电流幅值
    lineCurrentPrintToFile=Array{Float64}(undef,length(lines),points)
    println("!!!以下是各规划时刻各线路通过电流的真实值：")
    println(fileResult,"!!!以下是各规划时刻各线路通过电流的真实值：")
    for time in 1:points
        i=0
        for line in lines
            i+=1
            if checkLineAlive(line,lineData)
                #开方
                lineCurrentPrintToFile[i,time]=base_I*sqrt(abs(value(lineSquCurrent[(line[1],line[2],time)])))    
                println("****线路 $line 上流过的电流大小(A)为：",base_I*sqrt(abs(value(lineSquCurrent[(line[1],line[2],time)]))))
                println(fileResult,"****线路 $line 上流过的电流大小(A)为：",base_I*sqrt(abs(value(lineSquCurrent[(line[1],line[2],time)]))))
            else
                #ATTITION! 不可用的线路电流置为负无穷（这里用-9e9代替）
                lineCurrentPrintToFile[i,time]=-9e9
            end
        end
    end
    array2dPrintToFile(lineCurrentPrintToFile,lineCurrentFileResult)

    lineLossPrintToFile=Array{Float64}(undef,length(lines),points)
    println("!!!以下是各规划时刻网损真实值统计情况：")
    println(fileResult,"!!!以下是各规划时刻网损真实值统计情况：")
    for time in 1:points 
        i=0
        loss_T_System=0
        for line in lines
            i+=1
            if checkLineAlive(line,lineData)
                lineLossPrintToFile[i,time]=base_S/1000*rLineDict[line]*value(lineSquCurrent[(line[1],line[2],time)])
                println("****线路 $line 在 $time 规划时刻的网损大小(kW)为：",base_S/1000*rLineDict[line]*value(lineSquCurrent[(line[1],line[2],time)]))
                println(fileResult,"****线路 $line 在 $time 规划时刻的网损大小(kW)为：",base_S/1000*rLineDict[line]*value(lineSquCurrent[(line[1],line[2],time)]))
                loss_T_System+=base_S*rLineDict[line]/1000*value(lineSquCurrent[(line[1],line[2],time)])
            else
                #ATTITION!不可用线路网损置为负无穷（这里用-9e9代替）
                lineLossPrintToFile[i,time]=-9e9
            end
        end
        println("****在 $time 规划时刻的总体网损大小(kW)为：",loss_T_System)    
        println(fileResult,"****在 $time 规划时刻的总体网损大小(kW)为：",loss_T_System)    
    end
    array2dPrintToFile(lineLossPrintToFile,lineLossFileResult)

    if listMT!=()
        #ATTITION!行号仅代表listMT中元素的序号
        powerMTPrintToFile=Array{Complex}(undef,length(listMT),points)
        println("!!!以下是各MT在各规划时刻发出功率真实值统计情况：")
        println(fileResult,"!!!以下是各MT在各规划时刻发出功率真实值统计情况：")
        for time in 1:points
            i=0
            for mtNode in listMT
                i+=1
                if checkNodeFlagAlive(mtNode,nodeData,12)
                    powerMTPrintToFile[i,time]=complex(base_S/1000*(value(activePowerMT[(mtNode,time)]))
                    ,base_S/1000*(value(reactivePowerMT[(mtNode,time)])))
                    println("****MT $mtNode 在 $time 规划时刻的发出有功大小(kW)为：",base_S/1000*(value(activePowerMT[(mtNode,time)])))
                    println("****MT $mtNode 在 $time 规划时刻的运行成本为：",mtInfDict[mtNode][5]*(value(activePowerMT[(mtNode,time)])))
                    println(fileResult,"****含MT的节点 $mtNode 在 $time 规划时刻的发出有功大小(kW)为：",base_S/1000*(value(activePowerMT[(mtNode,time)])))
                    println(fileResult,"****MT $mtNode 在 $time 规划时刻的运行成本为：",mtInfDict[mtNode][5]*(value(activePowerMT[(mtNode,time)])))
                    println("****MT $mtNode 在 $time 规划时刻的发出无功大小(kVar)为：",base_S/1000*(value(reactivePowerMT[(mtNode,time)])))
                    println(fileResult,"****含MT的节点 $mtNode 在 $time 规划时刻的发出无功大小(kVar)为：",base_S/1000*(value(reactivePowerMT[(mtNode,time)])))
                else
                    #ATTITION!不可用MT出力置为0
                    powerMTPrintToFile[i,time]=0+0im
                end
            end 
        end
        array2dPrintToFile(powerMTPrintToFile,powerMTFileResult)
    end

    powerSubstationPrintToFile=Array{Complex}(undef,length(listSub),points)
    #ATTITION!行号仅代表listSub中元素的序号
    println("!!!以下是变电站主变出线功率情况")
    println(fileResult,"!!!以下是变电站主变出线功率情况")
    for time in 1:points
        i=0
        for subNode in listSub
            i+=1
            powerSubstationPrintToFile[i,time]=Complex(base_S/1000*(value(injectionActivePower[(subNode,time)])),base_S/1000*(value(injectionReactivePower[(subNode,time)])))
            println("****变电站 $subNode 在 $time 规划时刻从上级电网输入的有功大小(kW)为：",base_S/1000*(value(injectionActivePower[(subNode,time)])))
            println(fileResult,"****变电站 $subNode 在 $time 规划时刻从上级电网输入的有功大小(kW)为：",base_S/1000*(value(injectionActivePower[(subNode,time)])))
            println("****变电站 $subNode 在 $time 规划时刻从上级电网输入的无功大小(kVar)为：",base_S/1000*(value(injectionReactivePower[(subNode,time)])))
            println(fileResult,"****变电站 $subNode 在 $time 规划时刻从上级电网输入的无功大小(kVar)为：",base_S/1000*(value(injectionReactivePower[(subNode,time)])))
        end
    end
    array2dPrintToFile(powerSubstationPrintToFile,powerSubstationFileResult)

    println("!!!####测试用：显示b的情况：")
    for pair in ijtjit1Pair 
        println("****辅助变量b在 $pair 上的值为：" ,value(bAuxiliary[pair]))
    end
    
    loadEnergrizedPrintToFile=Array{Int64}(undef,numNode,points)
    println("!!!以下是负荷供电状态")
    println(fileResult,"!!!以下是负荷停电损失情况")
    for time in 1:points
        loadEnergrizedPrintToFile[1,time]=-999
        for node in rootFreeNodes
            loadEnergrizedPrintToFile[node,time]=round(Int64,value(nodeEnergized[(node,time)]))
            if loadEnergrizedPrintToFile[node,time]==0
                println("****警告！ 配电网 $node 号节点负荷已停电！")
                println(fileResult, "****警告！ 配电网 $node 号节点负荷已停电！")
            end
            println("****配电网 $node 号节点负荷供电状态为： ",loadEnergrizedPrintToFile[node,time])
            println(fileResult, "****配电网 $node 号节点负荷供电状态为： ",loadEnergrizedPrintToFile[node,time])
        end    
    end
    array2dPrintToFile(loadEnergrizedPrintToFile,loadEnergrizedFileResult)

    println("!!!以下是负荷停电损失情况")
    println(fileResult,"!!!以下是负荷停电损失情况")
    for node in nodes
        for t in 1:points
            println("****在 $t 规划时刻 节点 $node 负荷停电损失为： ",failureLoss*consumeP[node,t]*(1-value(nodeEnergized[(node,t)])))
            println(fileResult,"****在 $t 规划时刻 节点 $node 负荷停电损失为： ",failureLoss*consumeP[node,t]*(1-value(nodeEnergized[(node,t)])))
        end
    end
    
    println("##################结束##########################")
    println(fileResult,bus33Reconfiguration)

    println("******#测试用！！！ 显示x的情况***")
    for t in points
        for ij in lines
            println("****在 $t 规划时刻 线路$ij 的x值为： ",round(Int64,value(lineEnergized[(ij[1],ij[2],t)])))
            println(fileResult,"****在 $t 规划时刻 线路$ij 的x值为： ",round(Int64,value(lineEnergized[(ij[1],ij[2],t)])))
        end
    end

    println("******#测试用！！！ 显示y的情况***")
    for t in points
        for i in nodes
            println("****在 $t 规划时刻 节点$i 的y值为： ",round(Int64,value(nodeEnergized[(i,t)])))
            println(fileResult,"****在 $t 规划时刻 节点$i 的y值为： ",round(Int64,value(nodeEnergized[(i,t)])))
        end
    end

    println("******#测试用！！！ 显示gamma的情况***")
    for t in points
        for ij in listIC
            println("****在 $t 规划时刻 线路$ij 的gamma值为： ",round(Int64,value(closedLineDeEnergized[(ij[1],ij[2],t)])))
            println(fileResult,"****在 $t 规划时刻 线路$ij 的gamma值为： ",round(Int64,value(lineEnergized[(ij[1],ij[2],t)])))
        end
    end

    close(fileResult)
    close(switchOperationFileResult)
    close(powerMTFileResult)
    close(nodeVoltageFileResult)
    close(lineCurrentFileResult)
    close(lineLossFileResult)
    close(powerSubstationFileResult)
end

#运行主函数
dpSolverReconfiguraiton33Bus()



###############################
