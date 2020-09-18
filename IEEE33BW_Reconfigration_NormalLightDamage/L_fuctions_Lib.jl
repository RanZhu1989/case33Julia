#NOTE  2020.9.9 v0.2 
#NOTE  Author：Ran Zhu @ School of Cyber Engineering and Science, SEU
#NOTE  弹性配电网最优潮流库文件


##   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\  
##  $$  __$$\ $$  __$$\ $$  __$$\ $$  __$$\ $$ ___$$\ 
##  $$ /  \__|$$ /  \__|$$ /  \__|$$ /  \__|\_/   $$ |
##  $$ |      $$ |      $$ |      \$$$$$$\    $$$$$ / 
##  $$ |      $$ |      $$ |       \____$$\   \___$$\ 
##  $$ |  $$\ $$ |  $$\ $$ |  $$\ $$\   $$ |$$\   $$ |
##  \$$$$$$  |\$$$$$$  |\$$$$$$  |\$$$$$$  |\$$$$$$  |
##   \______/  \______/  \______/  \______/  \______/ 
##                                                    
##                                                    
##   

#TODO 下一步规划：用CSV DataFrame 包建立内存数据库提高JuMP模型参数计算速度

using DelimitedFiles
using DataStructures 

#设置JuMP框架下所需所有参数
function paraInit(DataPath,startPoint,endPoint)
    #Setting Time Point
    global points=endPoint-startPoint+1

    #Depacking DataPath
    loadFilePath=DataPath[1]
    pvFilePath=DataPath[2]
    wtFilePath=DataPath[3]
    paraLinePath=DataPath[4]
    paraNodePath=DataPath[5]
    pvFLAGpath=DataPath[6]
    mtFLAGpath=DataPath[7]
    wtFLAGpath=DataPath[8]
    lineFLAGpath=DataPath[9]
    pricePath=DataPath[10]

    #读取线路参数文件
    #reference:  M.E.Baran,1989,TPS
    readLineParaData=readdlm(paraLinePath,',',header=true)
    #每次调用求解器都会保存最后一步的动作作为下一次调用时开关的初始状态 故用global型 后同 
    global lineData=readLineParaData[1]
    global listLine=round.(Int64,lineData[:,1])
    global startNode=round.(Int64,lineData[:,2])
    global endNode=round.(Int64,lineData[:,3])
 
    #Create Parameters Array of Termination Lines
    global lineResistance=(lineData[:,4])/base_Z
    global lineImpedance=(lineData[:,5])/base_Z
    global lineSquCurrentMAX=(lineData[:,6])/base_I
    #读取联络开关位置，备用
    #global listFlagTieSwitch=lineData[:,7]
    
    #读取线路各时刻存活状态文件 [t,linenum]
    #线路存活时alpha才具有{0，1}属性，线路故障时强制断开，其alpha值为0
    readLineFLAGData=readdlm(lineFLAGpath,',',header=true)
    global lineFLAG=readLineFLAGData[1]

    #生成线路相关数据结构
    #总节点个数
    global numNode=max(maximum(startNode),maximum(endNode))
    #总线路条数
    global numLine=length(listLine)
    #遍历系统中所有节点的迭代器
    global nodes=1:numNode
    #遍历非根节点的迭代器
    #ATTITION!：这里人工指定了root=1
    global rootFreeNodes=2:numNode
    #(i,t) 节点-时间 二元组 t指规划尺度上的相对时刻 #ATTITION!在重载函数时需要重新设置
    global itPair=time1itPair(numNode,points)
    #生成标准情况下的线路(i,j)元组
    global lines=Tuple((startNode[i],endNode[i]) for i in 1:numLine)
    #生成(i,j)->NO.line的字典
    global lineDic=Dict(lines .=>(i for i in 1:numLine))
    #生成规划尺度上t时刻存活的线路(i,j,t)元组 #ATTITION!在重载函数时需要重新设置
    global linesAlive=findTupleAliveLinesTS(lines,lineData,lineFLAG,startPoint,endPoint)
    #生成规划尺度上t时刻故障线路(i,j,t)元组 #ATTITION!在重载函数时需要重新设置
    global linesFault=findTupleFaultLinesTS(lines,lineData,lineFLAG,startPoint,endPoint)
    #线路电阻值与阻抗的字典
    global rLineDict=Dict(lines .=>lineResistance)
    global xLineDict=Dict(lines .=>lineImpedance)
    #线路最大载流量的字典
    global maxCurrentLineDict=Dict(lines .=>lineSquCurrentMAX)
    #所有的(i,j,t)三元组穷举 t指规划时刻 #ATTITION!在重载函数时需要重新设置
    global time1NodePair=lines2TimePair(lines,points,1)
    global time0NodePair=lines2TimePair(lines,points,0)
    #各时刻存活的 (i,j,t) (j,i,t)三元组 t=1：points #ATTITION!在重载函数时需要重新设置
    global ijtjit1Pair=ijjiTime1Pair(linesAlive)
    #读取各线路初始状态 该初始状态stateInit将被solver.jl不断修改 
    global stateInit=round.(Int64,lineData[:,9])
    #ATTITION!在重载函数时需要重新设置
    #开关状态初始点，用于计算开关动作次数
    global lineInitStartPoint=ijPairInit(startNode,endNode,stateInit,numLine)
    #ATTITION!以下配置需要在Solver结尾处更新以供重载
    #线路通断start_value 用于给优化器设定初值增加收敛速度
    global alphaStartValue=Array{Int64}(undef,horizon-1,numLine)
    #注入有功、无功的start_value
    global injectionPStartValue=Array{Float64}(undef,horizon-1,numNode)
    global injectionQStartValue=Array{Float64}(undef,horizon-1,numNode)
    #首端有功无功的start_value
    global apparentPStartValue=Array{Float64}(undef,horizon-1,numLine)
    global apparentQStartValue=Array{Float64}(undef,horizon-1,numLine)
    #线路电流幅值平方的start_value
    global IsquStartValue=Array{Float64}(undef,horizon-1,numLine)
    #节点电压幅值平方的start_value
    global VsquStartValue=Array{Float64}(undef,horizon-1,numNode)

    #读取存活的MT PV WT (i,t) #ATTITION!在重载函数时需要重新设置
    global tsMTalive=creatEquipITpair("MT",points,DataPath,"Tuple")
    global tsPValive=creatEquipITpair("PV",points,DataPath,"Tuple")
    global tsWTalive=creatEquipITpair("WT",points,DataPath,"Tuple")
    global dicTSmtAlive=creatEquipITpair("MT",points,DataPath,"Dict")
    global dicTSpvAlive=creatEquipITpair("PV",points,DataPath,"Dict")
    global dicTSwtAlive=creatEquipITpair("WT",points,DataPath,"Dict")
    
    #读取节点参数文件
    readNodeData=readdlm(paraNodePath,',',header=true)
    global nodeData=readNodeData[1]
    #读取常开和常闭开关所在线路
    global listIC,listIO=readtopologyData(nodeData)
    #找出所有变电站节点元组
    global listSub=findSub(nodeData)
    #找出所有含PV节点元组
    global listPV=findPV(nodeData)
    #找出所有含WT节点元组
    global listWT=findWT(nodeData)
    #创建MT机组信息的字典
    #字典说明：键值=节点编号；值=(1有功下界，2有功上界，3无功下界，4无功上界，5运行成本)
    global mtInfDict=findPQmodeMT(nodeData,base_S)
    #找出所有含MT节点元组
    global listMT=Tuple(Set(keys(mtInfDict)))
    #找出普通负载节点元组
    global listPureLoad=findPureLoad(nodeData,numNode)
    #生成各节点的(pv,wt,mt,t)安装状态的元组 #ATTITION!在重载函数时需要重新设置
    global allocationDG=findAllocationDG(tsMTalive,tsPValive,tsWTalive,numNode,points)

    #读取价格参数 绝对时间尺度
    readPriceData=readdlm(pricePath,',',header=true)
    global priceData=readPriceData[1]
    #网上电价时间序列 
    global stGridPrice=Tuple(priceData[:,1])
    #馈电电价时间序列 备用
    #global stSellPrice=priceData[:,2]

    #读取TESTload #ATTITION!在重载函数时需要重新设置
    #ATTITION!.csv文件填入真实值，读取时转换为标幺值
    #已测试
    TESTLoadData=readdlm(loadFilePath,',')
    global consumeP=zeros(points,numNode)
    global consumeQ=zeros(points,numNode)
    for i in 1:33
        t=0
        for point in startPoint:endPoint
            t+=1
            #表格中单位是kW
            #TODO 把这些[t,i]顺序换一下
            global consumeP[t,i]=real(string2Complex(TESTLoadData[point,i]))*1000/base_S
            global consumeQ[t,i]=imag(string2Complex(TESTLoadData[point,i]))*1000/base_S  
        end
    end

    #读取TESTpv #ATTITION!在重载函数时需要重新设置
    #已测试
    readTESTPVData=readdlm(pvFilePath,',',header=true)
    TESTPVData=readTESTPVData[1]
    #识别每列PV安装的位置
    global pvbusNum=Tuple(readTESTPVData[2])
    #为后续计算注入功率方便，没安装PV的节点功率按0处理，WT相同
    global pvActivePower=zeros(points,numNode)
    global pvReactivePower=zeros(points,numNode)
    for i in listPV
        t=0
        for point in startPoint:endPoint
            t+=1
            #表格中单位是kW
            global pvActivePower[t,i]=real(string2Complex(TESTPVData[point,findall(x->x=="$i",pvbusNum)[1]]))*1000/base_S
            global pvReactivePower[t,i]=imag(string2Complex(TESTPVData[point,findall(x->x=="$i",pvbusNum)[1]]))*1000/base_S
        end
    end

    #读取TESTwt #ATTITION!在重载函数时需要重新设置
    #已测试
    readTESTWTData=readdlm(wtFilePath,',',header=true)
    TESTWTData=readTESTWTData[1]
    #识别每列WT安装的位置
    global wtbusNum=Tuple(readTESTWTData[2])
    global wtActivePower=zeros(points,numNode)
    global wtReactivePower=zeros(points,numNode)
    for i in listWT
        t=0
        for point in startPoint:endPoint
            t+=1
            #表格中单位是kW
            global wtActivePower[t,i]=real(string2Complex(TESTWTData[point,findall(x->x=="$i",wtbusNum)[1]]))*1000/base_S
            global wtReactivePower[t,i]=imag(string2Complex(TESTWTData[point,findall(x->x=="$i",wtbusNum)[1]]))*1000/base_S
        end
    end


end

function paraRolling(startPoint,endPoint,DataPath)
    #
    #Setting Time Point
    global points=endPoint-startPoint+1
    #(i,t) 节点-时间 二元组 t指规划尺度上的相对时刻 #ATTITION!在重载函数时需要重新设置
    global itPair=time1itPair(numNode,points)
    #生成规划尺度上t时刻存活的线路(i,j,t)元组 #ATTITION!在重载函数时需要重新设置
    global linesAlive=findTupleAliveLinesTS(lines,lineData,lineFLAG,startPoint,endPoint)
    #生成规划尺度上t时刻故障线路(i,j,t)元组 #ATTITION!在重载函数时需要重新设置
    global linesFault=findTupleFaultLinesTS(lines,lineData,lineFLAG,startPoint,endPoint)
    #所有的(i,j,t)三元组穷举 t指规划时刻 #ATTITION!在重载函数时需要重新设置
    global time1NodePair=lines2TimePair(lines,points,1)
    global time0NodePair=lines2TimePair(lines,points,0)
    #各时刻存活的 (i,j,t) (j,i,t)三元组 t=1：points #ATTITION!在重载函数时需要重新设置
    global ijtjit1Pair=ijjiTime1Pair(linesAlive)
    global lineInitStartPoint=ijPairInit(startNode,endNode,stateInit,numLine)
    #读取存活的MT PV WT (i,t) #ATTITION!在重载函数时需要重新设置
    global tsMTalive=creatEquipITpair("MT",points,DataPath,"Tuple")
    global tsPValive=creatEquipITpair("PV",points,DataPath,"Tuple")
    global tsWTalive=creatEquipITpair("WT",points,DataPath,"Tuple")
    #生成各节点的(pv,wt,mt)安装状态的矩阵[] #ATTITION!在重载函数时需要重新设置
    global allocationDG=findAllocationDG(tsMTalive,tsPValive,tsWTalive,numNode,points)

    #读取TESTload #ATTITION!在重载函数时需要重新设置
    #ATTITION!.csv文件填入真实值，读取时转换为标幺值
    #已测试
    TESTLoadData=readdlm(loadFilePath,',')
    global consumeP=zeros(points,numNode)
    global consumeQ=zeros(points,numNode)
    for i in 1:33
        t=0
        for point in startPoint:endPoint
            t+=1
            #表格中单位是kW
            #TODO 把这些[t,i]顺序换一下
            global consumeP[t,i]=real(string2Complex(TESTLoadData[point,i]))*1000/base_S
            global consumeQ[t,i]=imag(string2Complex(TESTLoadData[point,i]))*1000/base_S  
        end
    end

    #读取TESTpv #ATTITION!在重载函数时需要重新设置
    #已测试
    readTESTPVData=readdlm(pvFilePath,',',header=true)
    TESTPVData=readTESTPVData[1]
    #识别每列PV安装的位置
    global pvbusNum=Tuple(readTESTPVData[2])
    #为后续计算注入功率方便，没安装PV的节点功率按0处理，WT相同
    global pvActivePower=zeros(points,numNode)
    global pvReactivePower=zeros(points,numNode)
    for i in listPV
        t=0
        for point in startPoint:endPoint
            t+=1
            #表格中单位是kW
            global pvActivePower[t,i]=real(string2Complex(TESTPVData[point,findall(x->x=="$i",pvbusNum)[1]]))*1000/base_S
            global pvReactivePower[t,i]=imag(string2Complex(TESTPVData[point,findall(x->x=="$i",pvbusNum)[1]]))*1000/base_S
        end
    end

    #读取TESTwt #ATTITION!在重载函数时需要重新设置
    #已测试
    readTESTWTData=readdlm(wtFilePath,',',header=true)
    TESTWTData=readTESTWTData[1]
    #识别每列WT安装的位置
    global wtbusNum=Tuple(readTESTWTData[2])
    global wtActivePower=zeros(points,numNode)
    global wtReactivePower=zeros(points,numNode)
    for i in listWT
        t=0
        for point in startPoint:endPoint
            t+=1
            #表格中单位是kW
            global wtActivePower[t,i]=real(string2Complex(TESTWTData[point,findall(x->x=="$i",wtbusNum)[1]]))*1000/base_S
            global wtReactivePower[t,i]=imag(string2Complex(TESTWTData[point,findall(x->x=="$i",wtbusNum)[1]]))*1000/base_S
        end
    end
end




#############################################################################

function ijPairInit(startNode,endNode,stateInit,numLine)
    result=[]
    for i in 1:numLine
        push!(result,(startNode[i],endNode[i],stateInit[i]))
    end
    return  Tuple(result)
end

function string2Complex(complexString)
    #读取一个含a+b*im格式的字符串，然后转换为Complex类型
    indexPlus=findfirst(isequal('+'),complexString)
    indexIm=findfirst(isequal('i'),complexString)
    strReal=complexString[1:indexPlus-1]
    strImag=complexString[indexPlus+1:indexIm-1]
    return Complex(parse(Float64,strReal),parse(Float64,strImag))
end

function readtopologyData(titleFreeNodeDataArray)
    #输入去掉表头的NodeData
    #读取初始拓扑配置
    #返回常闭开关元组，常开开关元组
    listIC=[]
    listIO=[]
    #常闭
    indexIC=findall(x->x==1,titleFreeNodeDataArray[:,9])
    #常开
    indexIO=findall(x->x==0,titleFreeNodeDataArray[:,9])
    for i in 1:length(indexIC)
        push!(listIC,(round(Int64,titleFreeNodeDataArray[indexIC[i],2]),round(Int64,titleFreeNodeDataArray[indexIC[i],3])))
    end
    for i in 1:length(indexIO)
        push!(listIO,(round(Int64,titleFreeNodeDataArray[indexIO[i],2]),round(Int64,titleFreeNodeDataArray[indexIO[i],3])))
    end
    return Tuple(listIC), Tuple(listIO)
end


function array2dPrint2File(printArray,fileOpened)
    #将二维数组打印到文件
    #第一个参数是一个m*n矩阵 第二个参数是File=open(PATH,"w")形式
    for i in 1:(size(printArray))[1]
        for j in 1:(size(printArray))[2]
            print(fileOpened,printArray[i,j])
            print(fileOpened," ")
        end
        println(fileOpened)
    end
end

function checkNodeFlagAlive(i,NodeData,j)
    #已测试
    #检查一个节点是否存活 返回布尔量
    #存活指负荷\发电机是否退出
    #Data是node.csv文件，j=10 负载，j=11 PV j=12 MT
    if NodeData[i,j]==1
        return true
    else
        return false
    end
end

function checkLineAlive(line,lineData,t,lineFLAG)
    #已测试
    #检查一条线路(i,j)在t时刻是否存活
    #ATTITION! t指绝对时刻
    for i in 1:length(lineData[:,1])
        if (lineData[i,2],lineData[i,3])==line
            if lineFLAG[t,i]==1
                return true
            end        
        end
    end
    return false
end

function findSub(nodeData)
    #已测试
    #寻找并返回变电站节点编号的元组
    find=[]
    i=0
    for j in nodeData[:,2]
        i+=1
        if round.(Int64,j)==1
           push!(find,i) 
        end
    end
    return Tuple(find)
end

function findPV(nodeData)
    #已测试
    #寻找并返回PV节点编号的元组
    find=[]
    i=0
    for j in nodeData[:,3]
        i+=1
        if round.(Int64,j)==1
           push!(find,i) 
        end
    end
    return Tuple(find)
end

function findWT(nodeData)
    #已测试
    #寻找并返回WT节点编号的元组
    find=[]
    i=0
    for j in nodeData[:,4]
        i+=1
        if round.(Int64,j)==1
           push!(find,i) 
        end
    end
    return Tuple(find)
end

function findPureLoad(nodeData,nodeCount)
    #已测试
    #第二个输入node总个数
    #寻找并返回纯负载节点编号的元组
    Sub=Set(findSub(nodeData))
    PV=Set(findPV(nodeData))
    MT=Set(keys(findPQmodeMT(nodeData,base_S)))
    #并集运算
    temp=union(Sub,PV)
    temp=union(temp,MT)
    ALL=Set(i for i in 1:nodeCount)
    #差集运算
    Load=setdiff(ALL,temp)
    return Tuple(Load)
end

function findPQmodeMT(nodeData,base_S)
    #已测试
    #寻找并返回:
    #生成MT编号、运行成本、MT出力上下界的字典
    #ATTITION!生成的字典将实际参数值转换为标幺值
    #NOTE 这里假设MT运行在PQ模式，功率因数恒定，有容量限制
    findNode=[]
    findLP=[]
    findHP=[]
    findLQ=[]
    findHQ=[]
    findCost=[]
    i=0
    for j in nodeData[:,5]
        i+=1
        if round.(Int64,j)==1
           push!(findNode,i)
           #有功下界
           push!(findLP,(nodeData[i,6])/base_S)
           #有功上界 
           push!(findHP,(nodeData[i,7])/base_S) 
           #无功下界
           push!(findLQ,nodeData[i,8]/base_S)
           #无功上界
           push!(findHQ,(nodeData[i,9])/base_S) 
           #运行价格
           #ATTITION!在这里转换为标幺值下的价格
           push!(findCost,(nodeData[i,10])/(1000/base_S)) 
        end
    end
    TfindNode=Tuple(findNode)
    MTdata=Tuple((findLP[k],findHP[k],findLQ[k],findHQ[k],findCost[k]) for k in 1:length(TfindNode))
    MTdict=Dict(TfindNode .=>MTdata)
    return MTdict
end

function lines2TimePair(lines,points,start)
    #已测试
    #生成(i,j,t)组成的元组 其中(i,j)包含了所有线路 t是相对规划时刻 
    #ATTITION!  start==0 or start==1 ; points==endPoint-startPoint+1
    result=[]
    for line in lines
        for t in start:points
            #取元组中数据
            temp=(line[1],line[2],t)
            push!(result,temp)
        end
    end
    ijtPairTuple=Tuple(result)
    return ijtPairTuple
end

function ijjiTime1Pair(linesAlive)
    #已测试
    #将(i,j)组成的元组转换为(i,j,t),(j,i,t)组成的大型元组
    #将linesAlive中(i,j)颠倒后扩展 
    #ATTITION!  start==0 or start==1 ; points==endPoint-startPoint+1
    result=[]
    for line in linesAlive
        #取元组中数据
        temp1=(line[1],line[2],line[3])
        temp2=(line[2],line[1],line[3])
        push!(result,temp1)
        push!(result,temp2)
    end
    ijjitPairTuple=Tuple(result)
    return ijjitPairTuple
end

function time0ijtPair(lines)
    #已测试
    #生成(i,j,0)，
    result=[]
    for line in lines
        temp=(line[1],line[2],0)
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
    #ATTITION!  points==endPoint-startPoint+1
    temp=[i for i in 1:points]
    return Tuple(temp)
end

function time1itPair(numNode,points)
    #已测试
    #返回一个(i，t)组成的元组
    #ATTITION!  points==endPoint-startPoint+1
    result=[]
    for i in 1:numNode
        for t in 1:points
            temp=(i,t)
            push!(result,temp)
        end
    end
    return Tuple(result)
end

function findijNeighborNode(i,linesAlive,t)
    #已测试
    #返回节点i在t时刻与所有可用邻居组成的(i,j)元组
    #这里的可用指：相连的线路存活为真
    #i=节点标号 linesAlive=(i,j,t)元组集合 t=所检查的时刻
    result=[]
    for ijt in linesAlive
        if  t==ijt[3]
            if i in (ijt[1],ijt[2])
                push!(result,(ijt[1],ijt[2],t))
            end
        end
    end
    return Tuple(result)
end

function findIJJITNeighborNode(i,linesAlive)
    #已测试
    #ATTITION!这是针对(i,j)U(j,i)的一个无向搜索
    #ATTITION!  points==endPoint-startPoint+1
    #NOTE 某些情况下可能存在没有邻居的情况
    #NOTE 在重构问题中 都会存在至少一个邻居 （不形成孤岛）
    #返回节点i与所有可用邻居组成的(i,j,t)元组
    #这里的可用指：相连的线路存活为真
    #i=节点标号 lines=(i,j)元组集合 Data=读取的线路数据
    result=[]
    for ijt in linesAlive
        if i==ijt[1] 
            temp=(ijt[1],ijt[2],ijt[3])
            push!(result,temp)
        end
        if i==ijt[2]
            temp=(ijt[2],ijt[1],ijt[3])
            push!(result,temp)
        end
    end
    if result!=[]
        return Tuple(result)
    else
        return ()
    end
end

function findjkForwardNeighborPair(j,linesAlive,t)
    #已测试
    #返回节点j在t时刻给定方向与下游可形成连接的邻居节点组成的(j,k)元组
    #方向实际由csv文件的起讫点给出
    result=[]
    allPair=findijNeighborNode(j,linesAlive,t)
    for pair in allPair
        if pair[1]==j
            push!(result,pair)
        end
    end
    return Tuple(result)
end

function findijBackforwardNeighborPair(j,linesAlive,t)
    #已测试
    #返回节点j在给定方向与上游可形成连接的邻居节点组成的(i,j)元组
    #方向实际由csv文件的起讫点给出
    result=[]
    allPair=findijNeighborNode(j,linesAlive,t)
    for pair in allPair
        if pair[2]==j
            push!(result,pair)
        end
    end
    return Tuple(result)
end

function findTupleAliveLinesTS(lines,lineData,lineFLAG,startPoint,endPoint)
    #已测试
    #返回各时刻存活的线路组成的(i,j,t)元组
    result=[]
    for line in lines
        point=0
        for t in startPoint:endPoint
            point+=1
            if  checkLineAlive(line,lineData,t,lineFLAG)
                push!(result,(line[1],line[2],point))
            end
        end
    end
    return  Tuple(result)
end

function findTupleFaultLinesTS(lines,lineData,lineFLAG,startPoint,endPoint)
    #已测试
    #返回各故障的线路组成的(i,j,t)元组
    result=[]
    for line in lines
        point=0
        for t in startPoint:endPoint
            point+=1
            if  checkLineAlive(line,lineData,t,lineFLAG)==false
                push!(result,(line[1],line[2],point))
            end
        end
    end
    return  Tuple(result)
end

function nodeMTAliveitPair(listMT,dicMT,points)
    #返回存活的MT(node,t)
    #dicMT是MT字典
    #ATTITION!  points==endPoint-startPoint+1
    result=[]
    for node in listMT
        if dicMT[node][6]==1
            for t in 1:points
                push!(result,(node,t))
            end
        end
    end
    return  Tuple(result)
end

function findAllocationDG(tsMTalive,tsPValive,tsWTalive,numNode,points)
    #已测试
    #返回各时刻节点设备配置情况的数组
    result=Array{Tuple}(undef,points,numNode)
    # creat bin tuple (PV , WT ,MT ) for each （i,t)
    for t in 1:points  
        for node in 1:numNode
            result[t,node]=(1*in((node,t),tsPValive),1*in((node,t),tsWTalive),1*in((node,t),tsMTalive))
        end
    end
    return result
end

function creatEquipITpair(type,points,DataPath,mode)
    #已测试
    #创建特定类型设备的(i,t)组成的元组 即在t时刻i节点上存活
    #DataPath 6=PV 7=MT 8=WT
    if type=="PV"
        result=[]
        dictResult=[]
        readPVflagData=readdlm(DataPath[6],',',header=true)
        pvFLAG=readPVflagData[1]
        pvFLAGtitle=Tuple(readPVflagData[2])
        tempPVnode=[]
        for i in 1:length(pvFLAGtitle)
            push!(tempPVnode,parse(Int64,pvFLAGtitle[i]))
        end
        #创建一个列数->node号的字典
        num2NodeDic=SortedDict([i for i=1:length(pvFLAGtitle)].=>parse.(Int64,pvFLAGtitle))
        for t in 1:points
            for i in 1:length(pvFLAGtitle)
                if pvFLAG[t,i]==1
                    push!(result,(num2NodeDic[i],t))
                    push!(dictResult,(t,num2NodeDic[i]))
                end
            end
        end
        if mode=="Tuple"
            return Tuple(result)
        end
        if mode=="Dict"
            return Dict(Tuple(dictResult))
        end
    end
    if type=="MT"
        result=[]
        dictResult=[]
        readMTflagData=readdlm(DataPath[7],',',header=true)
        mtFLAG=readMTflagData[1]
        mtFLAGtitle=Tuple(readMTflagData[2])
        tempMTnode=[]
        for i in 1:length(mtFLAGtitle)
            push!(tempMTnode,parse(Int64,mtFLAGtitle[i]))
        end
        num2NodeDic=SortedDict([i for i=1:length(mtFLAGtitle)].=>parse.(Int64,mtFLAGtitle))
        for t in 1:points
            for i in 1:length(mtFLAGtitle)
                if mtFLAG[t,i]==1
                    push!(result,(num2NodeDic[i],t))
                    push!(dictResult,(t,num2NodeDic[i]))
                end
            end
        end
        if mode=="Tuple"
            return Tuple(result)
        end
        if mode=="Dict"
            return Dict(Tuple(dictResult))
        end 
    end
    if type=="WT"
        result=[]
        dictResult=[]
        readWTflagData=readdlm(DataPath[8],',',header=true)
        wtFLAG=readWTflagData[1]
        wtFLAGtitle=Tuple(readWTflagData[2])
        tempWTnode=[]
        num2NodeDic=SortedDict([i for i=1:length(wtFLAGtitle)].=>parse.(Int64,wtFLAGtitle))
        for i in 1:length(wtFLAGtitle)
            push!(tempWTnode,parse(Int64,wtFLAGtitle[i]))
        end
        for t in 1:points
            for i in 1:length(wtFLAGtitle)
                if wtFLAG[t,i]==1
                    push!(result,(num2NodeDic[i],t))
                    push!(dictResult,(t,num2NodeDic[i]))
                end
            end
        end
        if mode=="Tuple"
            return Tuple(result)
        end
        if mode=="Dict"
            return Dict(Tuple(dictResult))
        end   
    end
end

function mtNode2AliveDataColumn(node,DataPath)
    #返回mt节点对应与生存表中是第几列
    readMTflagData=readdlm(DataPath[7],',',header=true)
    mtFLAGtitle=Tuple(readMTflagData[2])
    return findfirst(isequal("$node"),mtFLAGtitle)
end

function pvNode2AliveDataColumn(node,DataPath)
    #返回pv节点对应与生存表中是第几列
    readPVflagData=readdlm(DataPath[6],',',header=true)
    pvFLAGtitle=Tuple(readPVflagData[2])
    return findfirst(isequal("$node"),pvFLAGtitle)
end

function wtNode2AliveDataColumn(node,DataPath)
    #返回wt节点对应与生存表中是第几列
    readWTflagData=readdlm(DataPath[8],',',header=true)
    wtFLAGtitle=Tuple(readWTflagData[2])
    return findfirst(isequal("$node"),wtFLAGtitle)
end

function findTimeFalutLine(t,linesFault)
    #返回t时刻所有的故障线路(i,j)
    result=[]
    for ijt in linesFault
        if  ijt[3]==t
            push!(result,(ijt[1],ijt[2]))
        end
    end
    return  Tuple(result)
end

function updateStateInit(switchOperationPrintToFile,points,stateInit,mode)
    #更新内存中的StateInit元组
    #mode="DP" 规划结果作为动作执行 initState更新为最后一个状态
    #mode="MPC" 仅规划结果的第一个动作执行 initState更新为第一个状态
    if mode=="DP"
        for i in size(switchOperationPrintToFile)[2]
            global stateInit[i]=switchOperationPrintToFile[points,i]
        end
    end
    if mode=="MPC"
        for i in size(switchOperationPrintToFile)[2]
            global stateInit[i]=switchOperationPrintToFile[1,i]
        end
    end
end