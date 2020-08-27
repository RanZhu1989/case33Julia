import Pkg
import MathOptInterface
using DelimitedFiles
using JuMP,MosekTools,Convex

####################静态数据集路径###########################
#ATTITION!现在填入的是测试用单步规划数据集
const loadFilePath="33BUSdata//forTEST//TEST_ts_33bus_load.csv"
const pvFilePath="33BUSdata//forTEST//TEST_ts_33bus_PV.csv"


####################ATTITION!以下为形成测试数据的语句，生成.jl文件时注释掉
#读取TESTload
#已测试
readTESTLoadData=readdlm(loadFilePath,',',header=true)
TESTLoadData=readTESTLoadData[1]
TESTlistP=TESTLoadData[:,3]
TESTlistQ=TESTLoadData[:,4]
global consumeP=zeros(33,9)
global consumeQ=zeros(33,9)
for iTESTLoad in 1:33
    for tTESTLoad in 1:points
        consumeP[iTESTLoad,tTESTLoad]=TESTlistP[iTESTLoad]
        consumeQ[iTESTLoad,tTESTLoad]=TESTlistQ[iTESTLoad]    
    end
end
#读取TESTPV
#已测试
readTESTPVData=readdlm(pvFilePath,',',header=true)
TESTPVData=readTESTPVData[1]
global PVP=zeros(33,9)
TESTlistPVP=TESTPVData[:,3]
for iTESTPV in 1:33
    for tTESTLoad in 1:points
        PVP[iTESTPV,tTESTLoad]=TESTlistPVP[iTESTPV]
    end
end



#######################全局常量##############################
#功率因数
const powerFactor=0.95
#规划时间尺度 26周
#points=4368
#单步规划 points=1
const points=1
#节点负荷时间序列矩阵 P or Q(i,t) 行i=节点号 列t=时间
#FIXME补全读取函数
#ATTITION!测试时不要运行下面这条语句
const consumeP,consumeQ=readConsumePower(loadFilePath,powerFactor,points)
#PV有功出力时间序列矩阵 P(i,t) 行i=节点号 列t=时间
#ATTITION!测试时不要运行下面这条语句
#FIXME补全读取函数
const PvP=readPVPower(pvFilePath,points)
#网损电价
const priceGrid=1#FIXME:网损成本需要设置，找参考文献！
#操作损耗
const costSwitch=1#FIXME:操作成本需要设置，找参考文献！
#节点电压上下界,原文中电压等级为12kV
const lowVNode=12000*0.9
const highVNode=12000*1.1
#每次开关最大动作次数
const maxNS=2
#无限制单步的优化用下面这个
#const maxNS=999
#大M 至少大于节点数即可
const bigM=999

#运行主函数
dpSolverReconfiguraiton33Bus()

###############################

function dpSolverReconfiguraiton33Bus()
    
    ####################读取必要数据文件########################################

    ##TODO 下面两个路径可作为参数传入
    #文件路径
    paraLinePath="33BUSdata//para_33bus_line.csv"
    paraNodePath="33BUSdata//para_33bus_node.csv" 

    #ATTITION!仅供测试用！形成脚本时注意注释掉
    paraLinePath="33BUSdata//forTEST//TEST_para_33bus_line.csv"
    paraNodePath="33BUSdata//forTEST//TEST_para_33bus_node.csv"

    #读取线路参数文件
    #reference:  M.E.Baran,1989,TPS
    readLineData=readdlm(paraLinePath,',',header=true)
    lineData=readLineData[1]
    listLine=round.(Int64,lineData[:,1])
    startNode=round.(Int64,lineData[:,2])
    endNode=round.(Int64,lineData[:,3])
    lineResistance=lineData[:,4]
    lineImpedance=lineData[:,5]
    lineCurrentMAX=lineData[:,6]
    #读取联络开关位置，备用
    listFlagTieSwitch=lineData[:,7]
    #生成线路相关数据结构
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
    #标准情况下逆序的线路(j,i)元组
    reverseLines=Tuple((endNode[i],startNode[i]) for i in 1:numLine)
    #线路电阻值与阻抗的字典
    rLineDict=Dict(lines .=>lineResistance)
    xLineDict=Dict(lines .=>lineImpedance)
    #线路最大载流量的字典
    maxCurrentLineDict=Dict(lines .=>lineCurrentMAX)
    #(i,j,t)三元组穷举，i,j∈line，t=0或1:points
    time0NodePair=linesToTimePair(lines,points,0)
    time1NodePair=linesToTimePair(lines,points,1)
    #(i,j,t) (j,i,t)三元组 t=1：points
    ijtjit1Pair=ijjiTimePair(lines,points,1)
    #(i,j,0)三元组
    ij0Pair=time0ijtPair(lines)
    #联络开关(i,j,0) 以及普通支路(i,j,0)
    ij0TiePair,ij0ComPair=ij0SET(lineData)
    
    #读取节点参数文件
    readNodeData=readdlm(paraNodePath,',',header=true)
    nodeData=readNodeData[1]
    #找出变电站节点元组
    listSub=findSub(nodeData)
    #找出PV节点元组
    listPV=findPV(nodeData)
    #找出并创建MT节点信息的字典
    #字典说明：键值=节点编号；值=(有功下界，有功上界，无功下界，无功上界，成本C)
    mtInfDict=findMT(nodeData)
    #找出MT节点元组
    listMT=Tuple(Set(keys(mtInfDict)))
    #找出普通负载节点元组
    listPureLoad=findPureLoad(nodeData,numNode)
    
    ###################建立JuMP模型##############################
    bus33Reconfiguration=Model(with_optimizer(Mosek.Optimizer))

    ###################目标函数######################################
    #ATTITION! 负载未定义 未编译
    @objective(bus33Reconfiguration,min,sum((sum((mtInfDict[mtNode][5])*(injectionActivePower[(mtNode,t)]
    +consumeP[mtNode,t]) for mtNode in listMT)
    +sum(priceGrid*sqrLineCurrent[(line[1],line[2],t)]*rLineDict[line] for line in lines)+
    sum(priceGrid*injectionActivePower[(subNode,t)] for subNode in listSub)
    +sum(costSwitch*numSwitchOperation[t])) for t in 1:points))
  
    ######################变量设置##################################
   
    #表示线路通断的alpha布尔矩阵
    #ATTITION!这里从时间0开始，代表初始状态
    @variable(bus33Reconfiguration,alpha[ijt in time0NodePair],Bin)

    #用于构建辐射状拓扑约束的中ST约束辅助变量b_{ijt}
    #ATTITION!这里从时间1开始
    #ATTITION!bAuxiliary包含(i,j)和它的反向对(j,i)
    @variable(bus33Reconfiguration,bAuxiliary[ijt in ijtjit1Pair],Bin)

    #用于构建辐射状拓扑约束的中SCF约束虚拟首端网络流变量apparentFictitiousFlow_{ijt}
    #ATTITION!这里从时间1开始
    @variable(bus33Reconfiguration,apparentFictitiousFlow[ijt in ijtjit1Pair])

    #用于构建辐射状拓扑约束的中SCF约束虚拟负载变量fictitiousLoad_{ijt}
    #ATTITION!这里从时间1开始
    @variable(bus33Reconfiguration,fictitiousLoad[it in itPair]==1)

    #用于构建开关次数约束的辅助变量lambda_{ijt}
    #ATTITION!这里从时间1开始
    @variable(bus33Reconfiguration,lambda[ijt in time1NodePair],Bin)

    #储存开关操作次数的辅助变量numSwitchOperation_{t}
    @variable(bus33Reconfiguration,numSwitchOperation[t in listPoints(points)]<=maxNS)

    #注入有功\无功功率 P Q _{i t}
   
    @variable(bus33Reconfiguration,injectionActivePower[it in itPair])
    @variable(bus33Reconfiguration,injectionReactivePower[it in itPair])

    #首端有功\无功功率 P Q _{i j t}
    @variable(bus33Reconfiguration,apparentActivePower[ijt in time1NodePair])
    @variable(bus33Reconfiguration,apparentReactivePower[ijt in time1NodePair])

    #l_{i j t}支路电流平方
    #每条支路最大载流量需要查lineData
    @variable(bus33Reconfiguration,sqrLineCurrent[ijt in time1NodePair])

    #节点电压幅值 V_{i t}
    @variable(bus33Reconfiguration,lowVNode<=nodeVoltage[it in itPair]<=highVNode)

    #######################补充的变量取值范围#########################

    #alpha的(i,j,0)值就是系统初始状态
    for pair in ij0TiePair
        @constraint(bus33Reconfiguration,alpha[pair]==0)
    end
    for pair in ij0ComPair
        @constraint(bus33Reconfiguration,alpha[pair]==1)
    end

    #每条支路最大载流量需要查lineData
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,sqrLineCurrent[(line[1],line[2],t)]<=maxCurrentLineDict[line])
        end
    end
   
    #注入功率：P_i=P发-P用
    #注意这里需要分以下情况：
    #1.纯负载节点或纯变电站节点;2.含PV的节点;3.含MT的节点
    #对于一般节点 
    #ATTITION!在PV Load数据读取函数编写完成前不要运行此部分！！
    for i in listPureLoad
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==-consumeP[i,t])
            @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==-consumeQ[i,t])
        end
    end       
    #对于含PV的节点
    for i in listPV
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==PvP[i,t]-consumeP[i,t])
            @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==-consumeQ[i,t])
    end
    #对于含MT的节点
    for i in listMT
        for t in 1:points
            mtInf=mtInfDict[i]
            @constraint(bus33Reconfiguration,mtInf[1]<=(injectionActivePower[(i,t)]+consumeP[(i,t)])<=mtInf[2])
            @constraint(bus33Reconfiguration,mtInf[3]<=(injectionReactivePower[(i,t)]+consumeQ[(i,t)])<=mtInf[4])
        end
    end

    ########约束条件########################
    #网络流模型的功率平衡方程
    #编译通过
    for node in nodes
        for t in 1:points
            if checkNodeAlive(node,nodeData)
                @constraint(bus33Reconfiguration,injectionActivePower[(node,t)]
                ==sum(apparentActivePower[(j,k,t)] for (j,k) in findjkForwardNeighborPair(node,lines,lineData) if j==node)-
                sum(apparentActivePower[(i,j,t)]-rLineDict[(i,j)]*sqrLineCurrent[(i,j,t)] 
                for (i,j) in findijBackforwardNeighborPair(node,lines,lineData) if j==node))
            
                @constraint(bus33Reconfiguration,injectionReactivePower[(node,t)]
                ==sum(apparentReactivePower[(j,k,t)] for (j,k) in findjkForwardNeighborPair(node,lines,lineData) if j==node)-
                sum(apparentReactivePower[(i,j,t)]-xLineDict[(i,j)]*sqrLineCurrent[(i,j,t)] 
                for (i,j) in findijBackforwardNeighborPair(node,lines,lineData) if j==node))
            end    
        end
    end

    #节点电压联系方程
    #编译通过
    for line in lines
        for t in 1:points
            if checkLineAlive(line,lineData)
                iVStart=line[1]
                jVEnd=line[2]
                @constraint(bus33Reconfiguration,(nodeVoltage[(iVStart,t)]-nodeVoltage[(jVEnd,t)])>=(-bigM*(1-alpha[(iVStart,jVEnd,t)])+
                2*(rLineDict[line]*apparentActivePower[(iVStart,jVEnd,t)]+xLineDict[line]*apparentReactivePower[(iVStart,jVEnd,t)]-
                ((rLineDict[line])^2+(xLineDict[line]^2))*sqrLineCurrent[(iVStart,jVEnd,t)])))

                @constraint(bus33Reconfiguration,(nodeVoltage[(iVStart,t)]-nodeVoltage[(jVEnd,t)])<=(bigM*(1-alpha[(iVStart,jVEnd,t)])+
                2*(rLineDict[line]*apparentActivePower[(iVStart,jVEnd,t)]+xLineDict[line]*apparentReactivePower[(iVStart,jVEnd,t)]-
                ((rLineDict[line])^2+(xLineDict[line]^2))*sqrLineCurrent[(iVStart,jVEnd,t)])))
            end
        end
    end

    #首端功率P^2+Q^2=VI^2的旋转二阶锥约束
    #TODO：JuMP提供的形式和标准有些不一样，需要进一步确认
    #编译通过
    for line in lines
        if checkLineAlive(line,lineData)
            iRSOCStart=line[1]
            jRSOCEnd=line[2]
            for t in 1:points
                @constraint(bus33Reconfiguration,[sqrLineCurrent[(iRSOCStart,jRSOCEnd,t)],0.5*nodeVoltage[(iRSOCStart,t)],
                apparentActivePower[(iRSOCStart,jRSOCEnd,t)],apparentReactivePower[(iRSOCStart,jRSOCEnd,t)]] in RotatedSecondOrderCone())
            end
        end
    end

    #基于ST+SCF生成树约束
    #Reference: ST = J.A.Taylor,2012,TPS; SCF = R.A.Jabr,2013,TPS
    #csv中的(start,end)中的start已被定义为end的父节点

    #ST约束
    #编译通过
    for line in lines
        if checkLineAlive(line,lineData)
            iSTStart=line[1]
            jSTEnd=line[2]
            for t in 1:points
                @constraint(bus33Reconfiguration,(bAuxiliary[(iSTStart,jSTEnd,t)]+bAuxiliary[(jSTEnd,iSTStart,t)])==alpha[(iSTStart,jSTEnd,t)])
            end
        end
    end
    for node in nodes
        for t in 1:points   
            #在寻找邻居节点时已自带线路存活检测功能
            @constraint(bus33Reconfiguration,sum(bAuxiliary[n] for n in findijtNeighborNode(node,lines,lineData,points))==1)
        end
    end
        #TODO：这里假设变电站就是唯一的root节点了
    subNeiijtPair=findijtNeighborNode(1,lines,lineData,points)
    for ijt in subNeiijtPair
        @constraint(bus33Reconfiguration,bAuxiliary[ijt]==0)
    end
    #SCF约束
    #编译通过
    for node in rootFreeNodes
        for t in 1:points
            if checkNodeAlive(node,nodeData)
                @constraint(bus33Reconfiguration,sum(apparentFictitiousFlow[(i,k,t)] for (i,k) in findjkForwardNeighborPair(node,lines,lineData) if i==node)+
                fictitiousLoad[(node,t)]==sum(apparentFictitiousFlow[(k,i,t)] for (k,i) in findijBackforwardNeighborPair(node,lines,lineData) if i==node))
            end
        end
    end
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,(apparentFictitiousFlow[(line[1],line[2],t)])<=(bigM*alpha[(line[1],line[2],t)]))
        end
    end
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,(apparentFictitiousFlow[(line[1],line[2],t)])>=(-bigM*alpha[(line[1],line[2],t)]))
        end
    end

end

##############################功能函数#############################################
function readConsumePower(filePath,powerFactor,points)
    #TODO：补全功能
    # 在输入时间尺度以及功率因数下返回各节点有功负载、无功负载.
    #格式：
    #
    return P,Q
end

function readPVPower(filePath,points)
    #TODO：补全功能
    # 在输入时间尺度以及功率因数下返回PV有功输出(假设该PV采用功率因数控制，不发无功).
    return P
end

function checkNodeAlive(i,Data)
    #已测试
    #检查一个节点是否存活 返回布尔量
    #存活指负荷\发电机是否退出
    #生存为假代表节点Pi=0
    if Data[i,9]==1
        return true
    else
        return false
    end
end

function checkLineAlive(line,Data)
    #已测试
    #检查一条线路(i,j)是否存活
    #存活为假说明该线路不得参与任何相关约束
    for i in 1:length(Data[:,1])
        if (Data[i,2],Data[i,3])==line
            if Data[i,8]==1
                return true
            else   
                return false
            end        
        end
    end
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
    MT=Set(keys(findMT(Data)))
    #并集运算
    temp=union(Sub,PV)
    temp=union(temp,MT)
    ALL=Set(i for i in 1:numNode)
    #差集运算
    Load=setdiff(ALL,temp)
    return Tuple(Load)
end

function findMT(Data)
    #已测试
    #寻找并返回:
    #生成MT编号、运行成本、MT出力上下界的字典
    findNode=[]
    findLP=[]
    findHP=[]
    findLQ=[]
    findHQ=[]
    findCost=[]
    i=0
    for j in Data[:,4]
        i+=1
        if round.(Int64,j)==1
           push!(findNode,i)
           push!(findLP,Data[i,5]) 
           push!(findHP,Data[i,6]) 
           push!(findLQ,Data[i,7]) 
           push!(findHQ,Data[i,8]) 
           push!(findCost,Data[i,9]) 
        end
    end
    TfindNode=Tuple(findNode)
    MTdata=Tuple((findLP[k],findHP[k],findLQ[k],findHQ[k],findCost[k]) for k in 1:length(TfindNode))
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
    #将(i,j)组成的元组转换为(i,j,t),(j,t,t)组成的大型元组
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

function findijtNeighborNode(i,lines,Data,points)
    #已测试
    #返回节点i与所有可用邻居组成的(i,j,t)元组
    #这里的可用指：相连的线路存活为真
    #i=节点标号 lines=(i,j)元组集合 Data=线路
    result=[]
    for line in lines
        if checkLineAlive(line,Data)
            if i==line[1]
                for t in 1:points
                    temp=(line[1],line[2],t)
                    push!(result,temp)
                end
            end
        end
    end
    return Tuple(result)
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

