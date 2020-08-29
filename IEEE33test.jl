#NOTE  2020.8.28 v0.5 
#NOTE  此.jl脚本用于测试锥规划模型下配网重构单步规划和长时间尺度动态规划数据生成器的性能
#NOTE  主要引用技术：
#NOTE  1. 基于网络流模型的Brach Flow Model : 1989, M.E.Baran, TPS
#NOTE  2. Brach Flow Model 相角松弛法 ：2013, M.Farivar, TPS
#NOTE  3. Spanning Tree 拓扑约束 : 2012, J.A.Taylor, TPS
#NOTE  4. Single Commodity Flow 拓扑约束： 2013, R.A.Jabr, TPS
#NOTE  5. 开关操作次数约束： 2016，Mohammad, TPS

import Pkg
import MathOptInterface
using DelimitedFiles
using JuMP,MosekTools,Convex

####################静态数据集路径###########################


#ATTITION!现在填入的是测试用单步规划数据集
#ATTITION!注意，原始数据需要在.jl脚本中转换为标幺值，后面计算采用标幺值计算

const loadFilePath="33BUSdata//forTEST//TEST_ts_33bus_load.csv"
const pvFilePath="33BUSdata//forTEST//TEST_ts_33bus_PV.csv"
const paraLinePath="33BUSdata//forTEST//TEST_para_33bus_line.csv"
const paraNodePath="33BUSdata//forTEST//TEST_para_33bus_node.csv"

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

#######################全局常量##############################
#基准容量
#用于计算标幺值
#Z_B=(U_B)^2 / S_B
#S_B = sqrt(3) * U_B * I_B
#r* =r/Z_B
#x* =x/Z_B
#p* =p/S_B
#q* =q/S_B

const base_S=100e6
const base_V=12e3
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
const priceLoss=0.5/(1000/base_S)#FIXME:网损成本需要设置，找参考文献！
#购电电价 单位monemy / kW*h 千瓦时
const priceGrid=0.6/(1000/base_S)#FIXME:购电电价需要设置，找参考文献！
#操作损耗
const costSwitch=0.5#FIXME:操作成本需要设置，找参考文献！
#节点电压上下界,原文中电压等级为12kV
#ATTITION!为线性化采用了电压幅值的平方作为变量！
const lowSquVNode=(0.9)^2
const highSquVNode=(1.1)^2
#ATTITION!每次开关最大动作次数 动态规划时使用这一条语句
#const maxNS=2
#ATTITION!无限制单步的优化用下面这条语句
#NOTE 设置为0 csv文件中修改上下限约束 可以作为潮流计算用
const maxNS=99
#大M 至少大于节点数即可
const bigM=33
####################ATTITION!以下为形成测试数据的语句，生成.jl文件时注释掉
#读取TESTload
#ATTITION!.csv文件填入真实值，读取时转换为标幺值
#已测试
readTESTLoadData=readdlm(loadFilePath,',',header=true)
TESTLoadData=readTESTLoadData[1]
TESTlistP=TESTLoadData[:,3]
TESTlistQ=TESTLoadData[:,4]
global consumeP=zeros(33,9)
global consumeQ=zeros(33,9)
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
global PVP=zeros(33,9)
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

function checkNodeAlive(i,Data)
    #已测试
    #检查一个节点是否存活 返回布尔量
    #存活指负荷\发电机是否退出
    #生存为假代表节点Pi=0
    if Data[i,10]==1
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
    i=0
    for j in Data[:,4]
        i+=1
        if round.(Int64,j)==1
           push!(findNode,i)
           push!(findLP,(Data[i,5])/base_S) 
           push!(findHP,(Data[i,6])/base_S) 
           push!(findLQ,(Data[i,7])/base_S) 
           push!(findHQ,(Data[i,8])/base_S) 
           push!(findCost,(Data[i,9])/(1000/base_S)) 
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
    lineResistance=(lineData[:,4])/base_Z
    lineImpedance=(lineData[:,5])/base_Z
    lineSquCurrentMAX=(lineData[:,6])/base_I
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
    #考虑线路生存状态的(i,j)U(j,i)
    ijjiLines=linesijTOlinesIJJI(lines,lineData)
    #线路电阻值与阻抗的字典
    rLineDict=Dict(lines .=>lineResistance)
    xLineDict=Dict(lines .=>lineImpedance)
    #线路最大载流量的字典
    maxCurrentLineDict=Dict(lines .=>lineSquCurrentMAX)
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
    mtInfDict=findMT(nodeData,base_S)
    #找出MT节点元组
    listMT=Tuple(Set(keys(mtInfDict)))
    #找出普通负载节点元组
    listPureLoad=findPureLoad(nodeData,numNode)
    
    ###################建立JuMP模型##############################
    bus33Reconfiguration=Model(with_optimizer(Mosek.Optimizer))
  
    ######################变量设置##################################
   
    #表示线路通断的alpha布尔矩阵
    #ATTITION!这里从时间0开始，代表初始状态
    @variable(bus33Reconfiguration,alpha[ijt in time0NodePair],Bin,base_name="线路通断Alpha")

    #用于构建辐射状拓扑约束的中ST约束辅助变量b_{ijt}
    #ATTITION!这里从时间1开始
    #ATTITION!bAuxiliary包含(i,j)和对应元素交换后的(j,i)
    @variable(bus33Reconfiguration,bAuxiliary[ijt in ijtjit1Pair],Bin,base_name="辅助变量b")

    #用于构建辐射状拓扑约束的中SCF约束虚拟首端网络流变量apparentFictitiousFlow_{ijt}
    #ATTITION!这里从时间1开始
    @variable(bus33Reconfiguration,apparentFictitiousFlow[ijt in ijtjit1Pair],base_name="虚拟网络流F")

    #用于构建辐射状拓扑约束的中SCF约束虚拟负载变量fictitiousLoad_{ijt}
    #ATTITION!这里从时间1开始
    @variable(bus33Reconfiguration,fictitiousLoad[it in itPair]==1,base_name="虚拟负载D")

    #用于构建开关次数约束的辅助变量lambda_{ijt}
    #ATTITION!这里从时间1开始
    @variable(bus33Reconfiguration,lambda[ijt in time1NodePair],Bin,base_name="开关变位辅助变量lambda")

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
    @variable(bus33Reconfiguration,lowSquVNode<=nodeSquVoltage[it in itPair]<=highSquVNode,base_name="节点电压幅值平方")

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
            @constraint(bus33Reconfiguration,lineSquCurrent[(line[1],line[2],t)]<=(maxCurrentLineDict[line])^2)
        end
    end
   
    #注入功率：P_i=P发-P用
    #注意这里需要分以下情况：
    #1.纯负载节点或纯变电站节点;2.含PV的节点;3.含MT的节点
    #对于一般节点 
    #编译通过
    for i in listPureLoad
        for t in 1:points
            @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==-consumeP[i,t])
            @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==-consumeQ[i,t])
        end
    end       
    #对于含PV的节点
    if listPV!=()
        for i in listPV
            for t in 1:points
                @constraint(bus33Reconfiguration,injectionActivePower[(i,t)]==PVP[i,t]-consumeP[i,t])
                @constraint(bus33Reconfiguration,injectionReactivePower[(i,t)]==-consumeQ[i,t])
            end
        end
    end
    
    #对于含MT的节点
    if listMT!=()
        for i in listMT
            for t in 1:points
                mtInf=mtInfDict[i]
                @constraint(bus33Reconfiguration,mtInf[1]<=(injectionActivePower[(i,t)]+consumeP[i,t])<=mtInf[2])
                @constraint(bus33Reconfiguration,mtInf[3]<=(injectionReactivePower[(i,t)]+consumeQ[i,t])<=mtInf[4])
            end
        end
    end
    

    ########约束条件########################
    #网络流模型的功率平衡方程
    #编译通过
    #已检查建模输出
    for node in nodes
        for t in 1:points
            if checkNodeAlive(node,nodeData)
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
    end

    #节点电压联系方程
    #编译通过
    for line in lines
        for t in 1:points
            if checkLineAlive(line,lineData)
                @constraint(bus33Reconfiguration,(nodeSquVoltage[(line[1],t)]-nodeSquVoltage[(line[2],t)])>=(-bigM*(1-alpha[(line[1],line[2],t)])
                +2*(rLineDict[line]*apparentActivePower[(line[1],line[2],t)]+xLineDict[line]*apparentReactivePower[(line[1],line[2],t)])
                -((rLineDict[line])^2+(xLineDict[line]^2))*lineSquCurrent[(line[1],line[2],t)]))

                @constraint(bus33Reconfiguration,(nodeSquVoltage[(line[1],t)]-nodeSquVoltage[(line[2],t)])<=(bigM*(1-alpha[(line[1],line[2],t)])
                +2*(rLineDict[line]*apparentActivePower[(line[1],line[2],t)]+xLineDict[line]*apparentReactivePower[(line[1],line[2],t)])
                -((rLineDict[line])^2+(xLineDict[line]^2))*lineSquCurrent[(line[1],line[2],t)]))
            end
        end
    end

    #首端功率P^2+Q^2=VI^2的旋转二阶锥约束
    #TODO：JuMP提供的形式和标准有些不一样，需要进一步确认
    #编译通过
    for line in lines
        if checkLineAlive(line,lineData)
            for t in 1:points
                @constraint(bus33Reconfiguration,[lineSquCurrent[(line[1],line[2],t)],0.5*nodeSquVoltage[(line[1],t)],
                apparentActivePower[(line[1],line[2],t)],apparentReactivePower[(line[1],line[2],t)]] in RotatedSecondOrderCone())
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
            for t in 1:points
                @constraint(bus33Reconfiguration,(bAuxiliary[(line[1],line[2],t)]+bAuxiliary[(line[2],line[1],t)])==alpha[(line[1],line[2],t)])
            end
        end
    end
    for node in rootFreeNodes
        for t in 1:points   
            #在寻找邻居节点时已自带线路存活检测功能
            @constraint(bus33Reconfiguration,sum(bAuxiliary[ijt] for ijt in rootFreeFindijtNeighborNode(node,ijjiLines,points))==1)
        end
    end
        #TODO：这里假设变电站就是唯一的root节点了
    subNeiijtPair=rootFreeFindijtNeighborNode(1,ijjiLines,points)
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

    #开关操作次数约束
    #Reference: Mohammad 2016, TPS
    #编译通过
    for line in lines
        for t in 1:points
            @constraint(bus33Reconfiguration,lambda[(line[1],line[2],t)]>=(alpha[(line[1],line[2],t)]-alpha[(line[1],line[2],t-1)]))
            @constraint(bus33Reconfiguration,lambda[(line[1],line[2],t)]>=(alpha[(line[1],line[2],t-1)]-alpha[(line[1],line[2],t)]))
        end
    end
    for node in nodes
        for t in 1:points
            @constraint(bus33Reconfiguration,sum(lambda[(i,j,t)] for (i,j) in findjkForwardNeighborPair(node,lines,lineData))<=maxNS)  
        end
    end
    
    ###################目标函数######################################
    #ATTITION! 负载未定义 
    #编译通过
    @objective(bus33Reconfiguration,Min,sum((sum((mtInfDict[mtNode][5])*(injectionActivePower[(mtNode,t)]
    +consumeP[mtNode,t]) for mtNode in listMT)
    +sum(priceLoss*lineSquCurrent[(line[1],line[2],t)]*rLineDict[line] for line in lines)
    +sum(priceGrid*injectionActivePower[(subNode,t)] for subNode in listSub)
    +sum(costSwitch*lambda[(line[1],line[2],t)] for line in lines)) for t in 1:points))
    
    
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

    println("!!!结束，优化最终状态标志位为：$state")
    println(fileResult,"!!!结束，优化最终状态标志位为：$state")
    println("!!!最优值:",objective_value(bus33Reconfiguration))
    println(fileResult,"!!!最优值:",objective_value(bus33Reconfiguration))

    #输出初始开关状态
    println("!!!初始各线路开关状态为: ")
    println(fileResult,"!!!初始各线路开关状态为: ")
    for line in lines
        if checkLineAlive(line,lineData)
            println("****线路 $line 的初始状态为：",Int64(value(alpha[(line[1],line[2],0)])))
            println(fileResult,"****线路 $line 的初始状态为：",Int64(value(alpha[(line[1],line[2],0)])))
        end
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
                switchOperationPrintToFile[n,time]=Int64(value(alpha[(line[1],line[2],time)]))
                println("****线路 $line 在 $time 规划时刻的状态为：",Int64(value(alpha[(line[1],line[2],time)])))
                println(fileResult,"****线路 $line 在 $time 规划时刻的状态为：",Int64(value(alpha[(line[1],line[2],time)])))
            else
                #ATTITION! 若线路当前不可用，状态置为-999
                switchOperationPrintToFile[n,time]=-999
            end
        end        
    end
    array2dPrintToFile(switchOperationPrintToFile,switchOperationFileResult)
    #println(switchOperationFileResult,switchOperationPrintToFile)
    println("!!!以下是各规划时刻发生动作的线路开关统计：")
    println(fileResult,"!!!以下是各规划时刻发生动作的线路开关统计：")
    for time in 1:points
        ns=0
        for line in lines
            if checkLineAlive(line,lineData)
                if value(lambda[(line[1],line[2],time)])==1
                    ns+=1
                    println("****线路 $line 上的断路器在 $time 规划时刻进行了1次操作")
                    println(fileResult,"****线路 $line 上的断路器在 $time 规划时刻进行了1次操作")
                end
            end
        end
        if ns==0
            println("****警告！配电网在 $time 规划时刻没有进行任何开关操作")
        end
    end

    #输出各节点电压幅值
    nodeVoltagePrintToFile=Array{Float64}(undef,length(nodes),points)
    println("!!!以下是各规划时刻各节点电压幅值的标幺值：")
    println(fileResult,"!!!以下是各规划时刻各节点电压幅值的标幺值：")
    for time in 1:points
        for node in nodes
            nodeVoltagePrintToFile[node,time]=sqrt(value(nodeSquVoltage[(node,time)]))
            println("****节点 $node 的电压幅值(pu.)为：",sqrt(value(nodeSquVoltage[(node,time)])))
            println(fileResult,"****节点 $node 的电压幅值(pu.)为：",sqrt(value(nodeSquVoltage[(node,time)])))
        end
    end
    array2dPrintToFile(nodeVoltagePrintToFile,nodeVoltageFileResult)
    #println(nodeVoltageFileResult,nodeVoltagePrintToFile)

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
                lineCurrentPrintToFile[i,time]=base_I*sqrt(value(lineSquCurrent[(line[1],line[2],time)]))    
                println("****线路 $line 上流过的电流大小(A)为：",base_I*sqrt(value(lineSquCurrent[(line[1],line[2],time)])))
                println(fileResult,"****线路 $line 上流过的电流大小(A)为：",base_I*sqrt(value(lineSquCurrent[(line[1],line[2],time)])))
            else
                #ATTITION! 不可用的线路电流置为0
                lineCurrentPrintToFile[i,time]=0
            end
        end
    end
    array2dPrintToFile(lineCurrentPrintToFile,lineCurrentFileResult)
    #println(lineCurrentFileResult,lineCurrentPrintToFile)

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
                #ATTITION!不可用线路网损置为0
                lineLossPrintToFile[i,time]=0
            end
        end
        println("****在 $time 规划时刻的总体网损大小(kW)为：",loss_T_System)    
        println(fileResult,"****在 $time 规划时刻的总体网损大小(kW)为：",loss_T_System)    
    end
    array2dPrintToFile(lineLossPrintToFile,lineLossFileResult)
    #println(lineLossFileResult,lineLossPrintToFile)

    if listMT!=()
        #ATTITION!行号仅代表listMT中元素的序号
        powerMTPrintToFile=Array{Complex}(undef,length(listMT),points)
        println("!!!以下是各MT在各规划时刻发出功率真实值统计情况：")
        println(fileResult,"!!!以下是各MT在各规划时刻发出功率真实值统计情况：")
        for time in 1:points
            i=0
            for mtNode in listMT
                i+=1
                if checkNodeAlive(mtNode,nodeData)
                    powerMTPrintToFile[i,time]=complex((base_S/1000*(value(injectionActivePower[(mtNode,time)])+consumeP[mtNode,time]))
                    ,(base_S/1000*(value(injectionReactivePower[(mtNode,time)])+consumeQ[mtNode,time])))
                    println("****MT $mtNode 在 $time 规划时刻的发出有功大小(kW)为：",base_S/1000*(value(injectionActivePower[(mtNode,time)])+consumeP[mtNode,time]))
                    println("****MT $mtNode 在 $time 规划时刻的运行成本为：",mtInfDict[mtNode][5]*base_S/1000*(value(injectionActivePower[(mtNode,time)])+consumeP[mtNode,time]))
                    println(fileResult,"****含MT的节点 $mtNode 在 $time 规划时刻的发出有功大小(kW)为：",base_S/1000*(value(injectionActivePower[(mtNode,time)])+consumeP[mtNode,time]))
                    println(fileResult,"****MT $mtNode 在 $time 规划时刻的运行成本为：",mtInfDict[mtNode][5]*base_S/1000*(value(injectionActivePower[(mtNode,time)])+consumeP[mtNode,time]))
                    println("****MT $mtNode 在 $time 规划时刻的发出无功大小(kVar)为：",base_S/1000*(value(injectionReactivePower[(mtNode,time)])+consumeQ[mtNode,time]))
                    println(fileResult,"****含MT的节点 $mtNode 在 $time 规划时刻的发出无功大小(kVar)为：",base_S/1000*(value(injectionReactivePower[(mtNode,time)])+consumeQ[mtNode,time]))
                else
                    #ATTITION!不可用MT出力置为0
                    powerMTPrintToFile[i,time]=0+0im
                end
            end 
        end
        array2dPrintToFile(powerMTPrintToFile,powerMTFileResult)
        #println(powerMTFileResult,powerMTPrintToFile)
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
    #println(powerSubstationFileResult,powerSubstationPrintToFile)

    # println("!!!####测试用：显示b的情况：")
    # for pair in ijtjit1Pair 
    #     println("****辅助变量b在 $pair 上的值为：" ,value(bAuxiliary[pair]))
    # end
    close(fileResult)
    close(switchOperationFileResult)
    close(powerMTFileResult)
    close(nodeVoltageFileResult)
    close(lineCurrentFileResult)
    close(lineLossFileResult)
    close(powerSubstationFileResult)

    println("##################结束##########################")
    #println(bus33Reconfiguration)
end

#运行主函数
dpSolverReconfiguraiton33Bus()



###############################
