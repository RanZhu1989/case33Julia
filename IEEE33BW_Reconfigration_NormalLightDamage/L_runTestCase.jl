#NOTE  2020.9.9 v0.9 
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
#Import OPFsolver: Adaptive or Normal Mode
include("L_opfSolverV2.jl")
#include("L_opfSolver_Normal.jl")
#Import Lib
include("L_fuctions_LibV2.jl")


#############################################################################################################
#Setting Data Path
loadFilePath=".//Xu2019TSG//Full_TEST//ISSDA_Load_FULL.csv"#1
pvFilePath=".//Xu2019TSG//Full_TEST//ts_PV_FULL.csv"#2
wtFilePath=".//Xu2019TSG//Full_TEST//ts_WT_FULL.csv"#3
paraLinePath=".//Xu2019TSG//Full_TEST//IEEE33BW_para_line.csv"#4
paraNodePath=".//Xu2019TSG//Full_TEST//IEEE38BW_para_node_Xu2019TSG.csv"#5
pvFLAGpath=".//Xu2019TSG//Full_TEST//FULL_FLAG_33bus_PV_Xu2019TSG.csv"#6
mtFLAGpath=".//Xu2019TSG//Full_TEST//FULL_FLAG_33bus_MT_Xu2019TSG.csv"#7
wtFLAGpath=".//Xu2019TSG//Full_TEST//FULL_FLAG_33bus_WT_Xu2019TSG.csv"#8
lineFLAGpath=".//Xu2019TSG//Full_TEST//FULL_FLAG_33bus_LINE_Xu2019TSG.csv"#9
pricePath=".//Xu2019TSG//Full_TEST//FULL_ts_33bus_Price_Xu2019TSG.csv"#10
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
global base_Z=(base_V)^2/(base_S)

#网损电价
global priceLoss=1/(1000/base_S)
#操作损耗
global costSwitch=10#FIXME:操作成本需要设置，找参考文献！
#停电单位损失
global failureLoss=10/(1000/base_S)#FIXME:停电成本需要设置，找参考文献！
global unitLossPenaltyCoefficient=1000#FIXME:停电成本需要设置，找参考文献！
#节点电压上下界,原文中电压等级为12kV
#ATTITION!为线性化采用了电压幅值的平方作为变量！
global lowSquVNode=(0.95)^2
global highSquVNode=(1.05)^2
#变电站主变出线电压幅值
global voltageSquSub=(1.00)^2
#假设MT具备黑启动能力，其节点电压可设置为1.00pu
global voltageBlackStartDG=1.00
#Minimum output of micro gas turbine
#zeroLowerBound=1 
#25% Pmax=0
ZLBflag=0
#joining in & dropping out control of MT PV WT
#MT JDC NOT available now
JDCflag=(0,0,0)
global maxNS=10
#大M 至少大于节点数即可
global bigM=35

#Setting MICOP optimizer parameters
#Setting relative gap termination tolerance
reGap=1e-2
#Setting maximum time spent by the mixed-integer optimizer, <0 => inf 
maxTime=200

#Setting planning start & end absolute time
#Total=12864 points
const planningStart=7
const planningEnd=7

#Setting absolute time pointer
#if using MPC or DP mode, horizon>=1
const horizon=1
global startPoint=planningStart
global endPoint=startPoint+horizon-1

#Mode="DP" or "MPC"
mode="MPC"
#request general report
global general_Report=1
############################################################################################
totalStartTime=now()
println("************Initialization.....StartingTime=$totalStartTime ************************")
#Init parameters for JuMP model
paraInit(DataPath,startPoint,endPoint)
#run optimizer for 1st time
println("************$(now())***** Initialization Fininsh！************************")
println("******************$(now())*****JuMP.jl->C++ Modeling************************************************")
dpSolverReconfiguraiton33Bus("Mosek","MPC",maxTime,reGap,JDCflag,ZLBflag)
#for test one-step mode
if planningStart==planningEnd
    totalEndTime=now()
    totalRunTime=parse(Float64,split(string(totalEndTime-totalStartTime)," ")[1])/1e3
    println("************Fininsh!*******************")
    println("************EndTime=$totalEndTime ************************")
    println("************Total Run Time= $totalRunTime seconds *******************")
else
    if horizon>=1
        #starting multi step planning
        for checkPoint in planningStart+1:planningEnd
            global startPoint=checkPoint
            if checkPoint+horizon-1<planningEnd
                global endPoint=startPoint+horizon-1
            else
                global endPoint=planningEnd
            end
            paraRolling(startPoint,endPoint,DataPath)
            dpSolverReconfiguraiton33Bus("Mosek","MPC",maxTime,reGap,JDCflag,ZLBflag)
        end
    end
    totalEndTime=now()
    totalRunTime=parse(Float64,split(string(totalEndTime-totalStartTime)," ")[1])/1e3
    println("************Fininsh!*******************")
    println("************EndTime=$totalEndTime ************************")
    println("************Total Run Time= $totalRunTime seconds *******************")

end
  

