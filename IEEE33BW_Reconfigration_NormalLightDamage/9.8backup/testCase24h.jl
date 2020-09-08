#NOTE  2020.9.6 v0.7 
#NOTE  Author：Ran Zhu @ School of Cyber Engineering and Science, SEU
#NOTE  此.jl文件用于评估正常运行及一般事故运行时自适应运行优化器的效果
#NOTE  **参数主要参考文献**：
#NOTE  1. IEEE33BW系统        1989, M.E.Baran, TPS
#NOTE  2. MT\PV\WT安装位置     2019, Xu, TSG
#NOTE  3. 负荷数据             Commission for Energy Regulation (CER). (2012). CER Smart Metering Project - Electricity Customer Behaviour Trial, 2009-2010 [dataset]. 1st Edition. Irish Social Science Data Archive. SN: 0012-00. www.ucd.ie/issda/CER-electricity
#NOTE  4. WT\PV地理及气象数据    https://www.renewables.ninja/ 
#NOTE  5. WT\PV数据生成模型     2016, Stefan, Energy  https://doi.org/10.1016/j.energy.2016.08.060.
#NOTE  6. 电缆参数             2011, Mancarella, IET-GT


# .___________. _______     _______.___________.     ______   .__   __.  __      ____    ____ 
# |           ||   ____|   /       |           |    /  __  \  |  \ |  | |  |     \   \  /   / 
# `---|  |----`|  |__     |   (----`---|  |----`   |  |  |  | |   \|  | |  |      \   \/   /  
#     |  |     |   __|     \   \       |  |        |  |  |  | |  . `  | |  |       \_    _/   
#     |  |     |  |____.----)   |      |  |        |  `--'  | |  |\   | |  `----.    |  |     
#     |__|     |_______|_______/       |__|         \______/  |__| \__| |_______|    |__| 

#Setting Current WorkPath
cd(@__DIR__)
#Import OPFsolver
include("opfSolver_LS9.6.jl")
#Import Lib
include("fuctions_Lib.jl")

#############################################################################################################
#Setting Data Path
loadFilePath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_ISSDA_536Days_TS_33BW_load.csv"#1
pvFilePath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_ts_33bus_PV_Xu2019TSG.csv"#2
wtFilePath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_ts_33bus_WT_Xu2019TSG.csv"#3
paraLinePath=".//33BUSdata//forTEST//Xu2019TSG//IEEE33BW_para_line.csv"#4
paraNodePath=".//33BUSdata//forTEST//Xu2019TSG//IEEE38BW_para_node_Xu2019TSG.csv"#5
pvFLAGpath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_FLAG_33bus_PV_Xu2019TSG.csv"#6
mtFLAGpath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_FLAG_33bus_MT_Xu2019TSG.csv"#7
wtFLAGpath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_FLAG_33bus_WT_Xu2019TSG.csv"#8
lineFLAGpath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_FLAG_33bus_LINE_Xu2019TSG.csv"#9
pricePath=".//33BUSdata//forTEST//Xu2019TSG//TEST24H_ts_33bus_Price_Xu2019TSG.csv"#10
#Packing Data Path to Tuple
#LOAD=1; pvSET=2; wtSET=3; paraLINE=4; paraNODE=5; pvFLAG=6; mtFLAG=7; wtFLAGpath=8; lineFLAG=9 price=10
DataPath=(loadFilePath,pvFilePath,wtFilePath,paraLinePath,paraNodePath,pvFLAGpath,mtFLAGpath,wtFLAGpath,lineFLAGpath,pricePath)
###################################################################################################################################

##############################################Setting Parameters###########################################################
#基准容量
#用于计算标幺值
#Z_B=(U_B)^2 / S_B
#S_B = sqrt(3) * U_B * I_B
#r* =r/Z_B
#x* =x/Z_B
#p* =p/S_B
#q* =q/S_B
#NOTE: 设置基准容量 => 设置为100MVA
global base_S=100e6
#NOTE: 设置基准电压 => 依据具体算例设置
global base_V=12.66e3
global base_I=(base_S)/sqrt(3)/base_V
global base_Z=(base_I)^2/base_S

#网损电价
global priceLoss=0.15
#操作损耗
global costSwitch=0.1#FIXME:操作成本需要设置，找参考文献！
#停电单位损失
global failureLoss=1/(1000/base_S)#FIXME:停电成本需要设置，找参考文献！
global unitLossPenaltyCoefficient=1000#FIXME:停电成本需要设置，找参考文献！
#节点电压上下界,原文中电压等级为12kV
#ATTITION!为线性化采用了电压幅值的平方作为变量！
global lowSquVNode=(0.95)^2
global highSquVNode=(1.05)^2
#变电站主变出线电压幅值
global voltageSquSub=(1.00)^2
#假设MT具备黑启动能力，其节点电压可设置为1.00pu
global voltageBlackStartDG=1.00
#ATTITION!每次开关最大动作次数 动态规划时使用这一条语句
#global maxNS=2
#ATTITION!无限制单步的优化用下面这条语句
#NOTE 设置为0 csv文件中修改上下限约束 可以作为潮流计算用
global maxNS=0
#大M 至少大于节点数即可
global bigM=35
#Setting DP time
global startPoint=13
global endPoint=14
############################################################################################
#generating parameters for JuMP model
paraInit(DataPath,startPoint,endPoint)
#run optimizer
dpSolverReconfiguraiton33Bus("Mosek")