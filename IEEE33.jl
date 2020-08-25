import Pkg
import MathOptInterface
using DelimitedFiles
using JuMP,MosekTools,Convex

function reconfiguraitonSolver()
    #读取必要数据文件########################################
    ##TODO 需要把该函数重新设计为所有参数可文件中读取，从而可形成自动化数据生成脚本
    #文件路径
    loadFilePath=""
    pvFilePath=""
    paraLinePath="33BUSdata//para_33bus_line.csv"
    paraNodePath="33BUSdata//para_33bus_node.csv" 
    
    #读取线路参数文件
    readLineData=readdlm(paraLinePath,',',header=true)
    lineData=readLineData[1]
    numLine=round.(Int64,lineData[:,1])
    startNode=round.(Int64,lineData[:,2])
    endNode=round.(Int64,lineData[:,3])
    lineResistance=lineData[:,4]
    lineImpedance=lineData[:,5]
    lineCurrentMAX=lineData[:,6]

    #生成线路相关数据结构
    #总节点个数，线路条数
    noNode=max(maximum(startNode),maximum(endNode))
    noLine=length(numLine)
    #遍历系统中所有节点的迭代器
    nodes=1:noNode
    #生成标准情况下的线路(i,j)元组
    lines=Tuple((startNode[i],endNode[i]) for i in 1:noLine)
    #线路电阻值与阻抗的字典
    rLineDict=Dict(lines .=>lineResistance)
    xLineDict=Dict(lines .=>lineImpedance)
    #线路最大载流量的字典
    maxCurrentLineDict=Dict(lines .=>lineCurrentMAX)
    
    #读取节点参数文件
    readNodeData=readdlm(paraNodePath,',',header=ture)
    nodeData=readLineData[1]
    #######################常量##############################
    #功率因数
    powerFactor=0.95
    #规划时间尺度 26周
    points=4368
    #微型燃气轮机单位功率的成本
    costMT=(6=>1,13=>1,22=>1,25=>1,28=>1,33=>1)#FIXME:MT的成本需要设置，找参考文献！
    #节点负荷时间序列
    consumeP,consumeQ=readConsumePower(loadFilePath,powerFactor,points)
    #PV有功出力时间序列
    PvP=readPVPower(pvFilePath,points)
    #网损电价
    priceGrid=1#FIXME:网损成本需要设置，找参考文献！
    #操作损耗
    costSwitch=1#FIXME:操作成本需要设置，找参考文献！
    #普通负荷节点集合
    setLoadNode=[2,3,4,5,8,9,10,11,12,14,15,16,17,19,20,21,23,24,26,27,29,30,31,32] 
    #含PV节点集合
    setPVNode=[7,18]
    #含MT节点集合
    setMTNode=[6,13,22,25,28,33]
    #变电站节点集合
    setSubstationNode=[1]
    #传输线参数 电阻real() ,阻抗imag()，最大传输电流
    #reference:  M.E.Baran,1989,TPS
    #FIXME:最大传输电流要重新设计，找参考文献！
    data
    parameterLine=Dict()
    parameterLine[(1,2)]=(0.0922+0.0470im,500)
    parameterLine[(2,3)]=(0.4930+0.2511im,500)
    parameterLine[(3,4)]=(0.3660+0.1864im,500)
    parameterLine[(4,5)]=(0.3811+0.1941im,500)
    parameterLine[(5,6)]=(0.8190+0.7070im,500)
    parameterLine[(6,7)]=(0.1872+0.6188im,500)
    parameterLine[(7,8)]=(0.7114+0.2351im,500)
    parameterLine[(8,9)]=(1.0300+0.7400im,500)
    parameterLine[(9,10)]=(1.0440+0.7400im,500)
    parameterLine[(10,11)]=(0.1966+0.0650im,500)
    parameterLine[(11,12)]=(0.3744+0.1238im,500)
    parameterLine[(12,13)]=(1.4680+1.1550im,500)
    parameterLine[(13,14)]=(0.5416+0.7129im,500)
    parameterLine[(14,15)]=(0.5910+0.5260im,500)
    parameterLine[(9,15)]=(2.0000+2.0000im,500)
    parameterLine[(15,16)]=(0.7463+0.5450im,500)
    parameterLine[(16,17)]=(1.2890+1.7210im,500)
    parameterLine[(17,18)]=(0.7320+0.5740im,500)
    parameterLine[(2,19)]=(0.1640+0.1565im,500)
    parameterLine[(19,20)]=(1.5042+1.3554im,500)
    parameterLine[(20,21)]=(0.4095+0.4784im,500)
    parameterLine[(8,21)]=(2.0000+2.0000im,500)
    parameterLine[(21,22)]=(0.7089+0.9373im,500)
    parameterLine[(12,22)]=(2.0000+2.0000im,500)
    parameterLine[(3,23)]=(0.4512+0.3083im,500)
    parameterLine[(23,24)]=(0.8980+0.7091im,500)
    parameterLine[(24,25)]=(0.8960+0.7011im,500)
    parameterLine[(25,29)]=(0.5000+0.5000im,500)
    parameterLine[(6,26)]=(0.2030+0.1034im,500)
    parameterLine[(26,27)]=(0.2842+0.1447im,500)
    parameterLine[(27,28)]=(1.0590+0.9337im,500)
    parameterLine[(28,29)]=(0.8042+0.7006im,500)
    parameterLine[(29,30)]=(0.5075+0.2585im,500)
    parameterLine[(30,31)]=(0.9744+0.9630im,500)
    parameterLine[(31,32)]=(0.3105+0.3619im,500)
    parameterLine[(32,33)]=(0.3410+0.5302im,500)
    parameterLine[(18,33)]=(0.5000+0.5000im,500)
    #节点电压上下界,原文中电压等级为12kV
    lowVNode=12000*0.9
    highVNode=12000*1.1
    #MT出力上下界
    lowMTactivePower=(6=>0,13=>0,22=>10,25=>0,28=>0,33=>0)
    highMTactivePower=(6=>99,13=>99,22=>99,25=>99,28=>99,33=>99)#FIXME:每台MT有功出力上界
    lowMTreactivePower=(6=>0,13=>0,22=>10,25=>0,28=>0,33=>0)
    highMTreactivePower=(6=>99,13=>99,22=>99,25=>99,28=>99,33=>99)#FIXME:每台MT无功出力上界
    #拓扑的时间序列
    originTopology=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),
    (10,11),(11,12),(12,13),(13,14),(14,15),(9,15),(15,16),(16,17),(17,18),(2,19),(19,20)
    ,(20,21),(8,21),(21,22),(12,22),(3,23),(23,24),(24,25),(25,29),(6,26),(26,27),
    (27,28),(28,29),(29,30),(30,31),(31,32),(32,33),(18,33)]
    networkTopology=zeros(33,33,points)
    networkTopology[1,2,1]=1
    networkTopology[2,3,1]=1
    networkTopology[3,4,1]=1
    networkTopology[4,5,1]=1
    networkTopology[5,6,1]=1
    networkTopology[6,7,1]=1
    networkTopology[7,8,1]=1
    networkTopology[8,9,1]=1
    networkTopology[9,10,1]=1
    networkTopology[10,11,1]=1
    networkTopology[11,12,1]=1
    networkTopology[12,13,1]=1
    networkTopology[13,14,1]=1
    networkTopology[14,15,1]=1
    networkTopology[9,15,1]=0
    networkTopology[15,16,1]=1
    networkTopology[16,17,1]=1
    networkTopology[17,18,1]=1
    networkTopology[2,19,1]=1
    networkTopology[19,20,1]=1
    networkTopology[20,21,1]=1
    networkTopology[8,21,1]=0
    networkTopology[21,22,1]=1
    networkTopology[12,22,1]=0
    networkTopology[3,23,1]=1
    networkTopology[23,24,1]=1
    networkTopology[24,25,1]=1
    networkTopology[25,29,1]=0
    networkTopology[6,26,1]=1
    networkTopology[26,27,1]=1
    networkTopology[27,28,1]=1
    networkTopology[28,29,1]=1
    networkTopology[29,30,1]=1
    networkTopology[30,31,1]=1
    networkTopology[31,32,1]=1
    networkTopology[32,33,1]=1
    networkTopology[18,33,1]=0
    #大M
    bigM=9999
 
    ########变量##################################
    #表示线路通断的alpha矩阵，t=0为初始态

   for oT in originTopology
       for t=1:points
            tempT=oT 
           
       end
   end

    bus33Reconfiguration=Model(with_optimizer(Mosek.Optimizer))
    @variable(bus33Reconfiguration,alpha[i123,j,t=1:points],Bin)
    #注入有功\无功功率
    injectionActivePower=zeros(33,points)
    @variable(bus33Reconfiguration,injectionActivePower,base_name="injP")
    injectionReactivePower=zeros(33,points)
    @variable(bus33Reconfiguration,injectionReactivePower,base_name="injQ")
    #首端有功\无功功率
    apparentActivePower=zeros(33,33,points)
    @variable(bus33Reconfiguration,apparentActivePower,base_name="apparentP")
    apparentReactivePower=zeros(33,33,points)
    @variable(bus33Reconfiguration,apparentReactivePower,base_name="apparentQ")
    #l_{ij}支路电流平方
    sqrLineCurrent=zeros(33,33,points)
    @variable(bus33Reconfiguration,sqrLineCurrent,base_name="sqrI")
    #节点电压幅值
    nodeVoltage=zeros(33,points)
    @variable(bus33Reconfiguration,nodeVoltage,base_name="nodeV")
    #辅助变量b_{ij}
    isParent=zeros(33,33,points)
    @variable(bus33Reconfiguration,isParent,base_name="b")
    #储存开关操作次数
    numSwitchOperation=zeros(points)
    @variable(bus33Reconfiguration,numSwitchOperation,base_name="NS")
    #动作次数约束所用到的辅助布尔型变量
    lambda=zeros(33,33,points)
    @variable(bus33Reconfiguration,lambda,base_name="lambda")

   
    ########变量取值范围#########################
    
    #布尔量
    for i in networkTopology
        @constraint(bus33Reconfiguration,i in MOI.ZeroOne())        
    end

    for i in isParent
        @constraint(bus33Reconfiguration,i in MOI.ZeroOne())
    end

    for i in lambda
        @constraint(bus33Reconfiguration,i in MOI.ZeroOne())      
    end
    #注入功率：
    #对于一般节点 
    for i in setLoadNode 
        for t=1:points
            @constraint(bus33Reconfiguration,injectionActivePower[i,t]==-consumeP[i,t])
            @constraint(bus33Reconfiguration,injectionReactivePower[i,t]==-consumeQ[i,t])
        end       
    end
    #对于含PV的节点
    for i in setPVNode
        for t=1:points
            @constraint(bus33Reconfiguration,injectionActivePower[i,t]==PvP[i,t]-consumeP[i,t])
            @constraint(bus33Reconfiguration,injectionReactivePower[i,t]==-consumeQ[i,t])
    end
    #对于含MT的节点
    for i in setMTNode
        for t=1:points
            @constraint(bus33Reconfiguration,lowMTactivePower<=(injectionActivePower[i,t]+consumeP[i,t])<=highMTactivePower)
            @constraint(bus33Reconfiguration,lowMTreactivePower<=(injectionReactivePower[i,t]+consumeQ[i,t])<=highMTreactivePower)
    end
    #传输线电流
    for k in keys(parameterLine)
        for t=1:points 
            l=parameterLine[k]
            @constraint(bus33Reconfiguration,sqrLineCurrent[i,j,t]<=l[2])  
        end        
    end
    #节点运行电压
    for i=1:33
        for t=1:points
            @constraint(bus33Reconfiguration,lowVNode<=nodeVoltage[i,t]<=highVNode)   
        end 
    end
    #每步规划时开关操作不多于两次
    for t=1:points
        @constraint(bus33Reconfiguration,numSwitchOperation[t]<=2)       
    end
    
 

    ########约束条件########################
    #网络流模型的功率平衡方程
    for t=1:points
        for node in checkNodes(networkTopology,t)
            tempP1=0
            tempQ1=0
            for pair in findNodeChild(node)
                tempP1+=apparentActivePower[pair]
                tempQ1+=apparentReactivePower[pair]
            end
            tempP2=0
            tempQ2=0
            for pair in findNodeParent(node)
                temp=parameterLine[pair]
                tempP2+=apparentActivePower[pair]-real(temp[1])*sqrLineCurrent[pair(1),pair(2),t]
                tempQ2+=apparentReactivePower[pair]-imag(temp[1])*sqrLineCurrent[pair(1),pair(2),t]
            end
            @constraint(bus33Reconfiguration,lowVNode<=nodeVoltage[i,t]<=highVNode)
            @constraint(bus33Reconfiguration,lowVNode<=nodeVoltage[i,t]<=highVNode)
        end
    end


    #目标函数
    @objective(bus33Reconfiguration,min,)








end

################辅助函数#################################
function readConsumePower(filePath,powerFactor,points)
    #TODO：补全功能
    # 在输入时间尺度以及功率因数下返回各节点有功负载、无功负载.
    return P,Q
end

function readPVPower(filePath,points)
    #TODO：补全功能
    # 在输入时间尺度以及功率因数下返回PV有功输出(假设该PV采用功率因数控制，不发无功).
    return P
end

function checkNodes(Link,t)
    #TODO：补全功能
    #输入拓扑及时刻，返回该拓扑下节点的列表
    return nodeList
end

function findNodeChild(j)
    #TODO：补全功能
    #输入一个节点编号，返回 自身与所有近邻的下游节点组成元组的列表
    return childTuple
end

function findNodeParent(i)
    #TODO：补全功能
    #输入一个节点编号，返回 自身与所有近邻的上游节点组成的元组列表
    return parentTuple
end