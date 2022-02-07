function Func_birthdeath_GC_logistic6(BindFlag,paramFlag,NAg,IC1,IC2,ClusterSize,A,CT,ProtF,OnTargetRCF,OffTargetRCF,FitAdv,ProtAg,id,memnumber,Corr12,Multiplier,Agcutoff)
rng(1000)
% AllClones = [1:1:228];
expon=1;
global expon;
AllClones = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,101,102,103,104,105,106,107,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,133,136,137,138,139,141,142,143,144,145,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,190,195,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228];

OnTargetClonesNPSS = [11 13 14 152 167 169 170 171 172 173]; %NPSS
AllClonesNPSS = [1 3 4 5 6 7 8 9 10 11 12 13 14 91 139 141 144 145 147 148 149 150 151 152 153 155 156 157 158 159 160 161 162 163,...
    164 165 166 167 168 169 170 171 172 173 174 175 190 195 209 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228];

OffTargetClonesNPSS = setdiff(AllClonesNPSS,OnTargetClonesNPSS);

OnTargetClonesNPHA = [11 13 14 128 129 152 153 167 169 170 171 172 173]; % NPHA

% OffTargetClonesNPHA = setdiff([1:1:228],AllClonesNPSS);
OffTargetClonesNPHA = setdiff(AllClones,AllClonesNPSS);


% The parameter is the antigen

% 1 : HA
% 2 : NPHA
% 3 : NPSS
% 4 : 40VirusHA

switch ProtAg
    case 1334
        NP1s = '1HA';NP2s = 'NPSS';NP3s = 'NPSS';NP4s = '40VirusHA';
    case 3334
        NP1s = 'NPSS';NP2s = 'NPSS';NP3s = 'NPSS';NP4s = '40VirusHA';
    case 2224
        NP1s = 'NPHA';NP2s = 'NPHA';NP3s = 'NPHA';NP4s = '40VirusHA';
    case 2334
        NP1s = 'NPHA';NP2s = 'NPSS';NP3s = 'NPSS';NP4s = '40VirusHA';
    case 4334
        NP1s = '40VirusHA';NP2s = 'NPSS';NP3s = 'NPSS';NP4s = '40VirusHA';
    case 4444
        NP1s = '40VirusHA';NP2s = '40VirusHA';NP3s = '40VirusHA';NP4s = '40VirusHA';
    case 1114
        NP1s = '1HA';NP2s = '1HA';NP3s = '1HA';NP4s = '40VirusHA';
	case 4344
	NP1s = '40VirusHA';NP2s = 'NPSS';NP3s = '40VirusHA';NP4s = '40VirusHA';

end

ProtAg
NP1s
NP2s
NP4s


q01Arm=2;
q02Arm=2;

AgThExp = 100;

ImmunClass = ImmunogenicityAgClass(100,100);

% NP = '1HA';
NP = NP1s;



KonArr = KonEpitopes(NP,q01Arm,q02Arm);

NAg_T2 = NAg;
NAg_T3 = NAg;
% CT a parameter describing the number of Tfh cells in a GC and thus the
% competitiveness of the process.
NAg = single(NAg);
NAg_T2 = single(NAg_T2);
NAg_T3 = single(NAg_T3);
rng('shuffle');
dt=0.1;
lambda = single(0.1);
lambda0 = single(1.5);
lambda0Exp = 0.42;
% lambda0Exp = 0;
mu = 0.6;
muExpDecay = 0.01;
% muExpDecay = 0;

Diff = single(0.05);
runNum = 14;
Nsp = KonArr.EpitopeNum
NumberOfCells4Seed = 200;
% NumberOfCells4Seed = 1;
T0 = 1;
T1 = 16;
T2 = 16;
T3 = 0;
T4 = 0;
% % 
% T0 = 1;
% T1 = 3;
% T2 = 3;
% T3 = 3;
% T4 = 0;
% 
% T0 = 1;
% T1 = 1;
% T2 = 1;
% T3 = 1;
% T4 = 0;

% TMemoryGrowth = 5;
TMemoryGrowth = T1+T0; %note I changed to 0 for the boost later.

% TExpansionDecay = 1;
% TExpansionDecay = (T1+T0)-TMemoryGrowth; % The decay of the memory B cells undergoing the expansion
% TExpansionDecayAgID = (T1+T0);% The decay of the memory B cells which do not react to the Ag

TExpansionDecay = 150; %normally 250
TExpansionDecayAgID = 0;

TMemoryGrowthFluchallenge = 0;

ICMemExpMultiplier = 4; % The B cells that expand are exposed to more Ag than B cells in the GC.

FCFlag = 0; % Flu Challange flag

Nmax = 5000;
FlagDeathStop = 0;
Jackpotflag = 0;
JPThreshold = 1.5;
BondsNum = 1;
NBCR = 20;
SecBind = 2;

qMFPTImmuneComp = 0.1;
q1=single([0.1 0.1 0.1 0.1]);

q2=single([0.05 0.05]);
tau = 1;
% tau = 5;
EAbAg = single([-0.5 0.1 0.1]);
EAgMem = single(2);
qMutP = 0.9;
% wcv = [0.05 0.95];
wcv = [0 1];
% wcv = [0.1 0.9]; % the relative contribution to the total affinity of the conserved resiude wcv(1) and the variable one wcv(2)
rcv = [0.1 0.9]; % the contribution to the affinity of the conserved part vs. the variable part

%%%%%%  Extract parameter %%
ExtcFlag = 1; % extract the whole immune complex
% ExtcFlag = 2; % gradual extraction

%%%%%%%% Remove B cells during the GCR %%%
% ExitProb = 0.05;
ExitProb = 0.1;
% ExitProb = 0.5;

EnergySeedClones = 0;
EnterProbAboveMed = 1;
EnterProbBelowMed = 0;
% AgThresholdGCSeed = 200;

NAg
% Pvc=0.5;
% Pvv=0.5;
OnTargetRCF
OffTargetRCF
FitAdv


CellsPerCloneArr = zeros(Nsp,1);
if(strcmp(NP1s,'NPSS'))
    OnTargetClones = ClonesList('NPSS');
    for n=1:Nsp
        if(find(OnTargetClones==n))
            CellsPerCloneArr(n) = OnTargetRCF;
        else
            CellsPerCloneArr(n) = OffTargetRCF;
        end
    end
    elseif(strcmp(NP1s,'NPHA') & strcmp(NP2s,'NPHA'))
    OnTargetClones = ClonesList('NPHA');
    for n=1:Nsp
        if(find(OnTargetClones==n))
            CellsPerCloneArr(n) = OnTargetRCF;
%             CellsPerCloneArr(n) = OnTargetRCF;
        elseif(find(OffTargetClonesNPSS==n))
%             CellsPerCloneArr(n) = OffTargetRCF;
            CellsPerCloneArr(n) = OffTargetRCF;
        elseif(find(OffTargetClonesNPHA==n))
%             CellsPerCloneArr(n) = 5;
            CellsPerCloneArr(n) = OffTargetRCF;
        end
        
    end
    elseif(strcmp(NP1s,'40VirusHA') & strcmp(NP2s,'40VirusHA'))
    OnTargetClones = ClonesList('NPHA');
    for n=1:Nsp
        if(find(OnTargetClones==n))
            CellsPerCloneArr(n) = OnTargetRCF;
        elseif(find(OffTargetClonesNPSS==n))
            CellsPerCloneArr(n) = OffTargetRCF;
        elseif(find(OffTargetClonesNPHA==n))
            CellsPerCloneArr(n) = OffTargetRCF;
        end
        
    end
    elseif(strcmp(NP1s,'1HA') & strcmp(NP2s,'1HA'))
    OnTargetClones = ClonesList('NPHA');
    for n=1:Nsp
        if(find(OnTargetClones==n))
            CellsPerCloneArr(n) = OnTargetRCF;
        elseif(find(OffTargetClonesNPSS==n))
            CellsPerCloneArr(n) = OffTargetRCF;
        elseif(find(OffTargetClonesNPHA==n))
            CellsPerCloneArr(n) = OffTargetRCF;
        end
        
    end
else
    OnTargetClones = ClonesList('NPSS');
    for n=1:Nsp
        if(find(OnTargetClones==n))
            CellsPerCloneArr(n) = OnTargetRCF;
        elseif(find(OffTargetClonesNPSS==n))
            CellsPerCloneArr(n) = OffTargetRCF;
        elseif(find(OffTargetClonesNPHA==n))
            CellsPerCloneArr(n) = OffTargetRCF;
        end
        
    end
    

    %%%%% Initial fraction of cells per clone %%%%%%
%     CellsPerCloneArr(OffTargetClonesNPHA) = OffTargetRCF;
%    CellsPerCloneArr(OffTargetClonesNPSS) = 1;
%    CellsPerCloneArr(OnTargetClonesNPSS) = 7;
   
% CellsPerCloneArr(OffTargetClonesNPSS) = 10;
% CellsPerCloneArr(OnTargetClonesNPSS) = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



% for n=1:Nsp
%     if(find(OnTargetClones==n))
%         CellsPerCloneArr(n) = OnTargetRCF;
%     else
%         CellsPerCloneArr(n) = OffTargetRCF;
%     end
% end

% CellsPerCloneArr = zeros(Nsp,1);
% for n=1:Nsp
%     if(find(OnTargetClones==n))
%         CellsPerCloneArr(n) = OnTargetRCF;
%     else
%         CellsPerCloneArr(n) = OffTargetRCF;
%     end
% end

%%%%%%%%%%%%
%%% magnitude of boost from interaction with HA
bOff = 1;
bOn = 1;
%%%%%%%%%%%
%%% T cells epitope specific CT boost
COff = 1
COn = 1
%%%
%%%%%%%%%%%


CorrVC = 0;
Corr23=Corr12;
Corr13=Corr12*Corr23 - sqrt(Corr12^2 * Corr23^2 -Corr12^2 - Corr23^2 +1);
CorrVV=[1 Corr12 Corr13; Corr12 1 Corr23; Corr13 Corr23 1];
MutVC = 0;
MutVV = 0;
NAnt=length(CorrVV);

ICTypeNums = [IC1 IC2]
% ProtF is the protocol flag
switch ProtF
    case 1
        AgICprop = [0.5 0.5 ; 0.5 0.5]
    case 2
        AgICprop = [1 0 ; 0.5 0.5]
    case 3
        AgICprop = [1 0 ; 0 1]
    case 4
        AgICprop = [1 0 ; 1 0]
end

% AgICprop = [0 1 ; 1 0]
% binding to different Ag molecules
% BindFlag=1 this is an Ag density simulation. only one kon
% BindFlag=2 this is an Ag density simulation. There are two kon
% BindFlag=5 this is an Ag density simulation. only one kon ; initial state is 00
% binding to the same Ag molecule
% BindFlag=3 This is an Ag density simulation. only one kon.
% BindFlag=4 This is an Ag density simulation. There are two kon. the first
% BindFlag=6 This is an Ag density simulation. only one kon ; initial state is 00
% is Ag density dependent and the second is not.

% BindFlag=8 This is an Ag density simulation. only one kon ; initial state
% is 00. There are two parameters: MAg and p (A) from which NAg is
% randomize.

% BindFlag = 4;

% paramFlag = 3;
% paramFlag=1 Mutate EAbAg
% paramFlag=2 Mutate q1
% paramFlag=3 Mutate q2
% paramFlag=4 Mutate q2

% ClusterSize = 5;

% when BindFlag=1 this is an Ag density simulation. when BindFlag=0 it's the double binding on the same Ag molecule % for 2 we consider q1 and q2 % for 3 there is cooperativity
a=1; % radius covered by the Ab
% A=30; % radius of the synapse
% NAg = 1; % number of Ag mulecules in radius A
F = 1e-12; % The force at which the B cell pulls the Ag-IC (immune complex)
% F = 0; % The force at which the B cell pulls the Ag-IC (immune complex)


%%%%%%%%% Initial density / k on dependency

% KonArr = KonDensity(FitAdv);

%%% Number of IC of each type

global nalpha nbeta
global ICArr;

clear AncArr;
clear bcArr;
clear cellnumGrowth;
clear cellnumSel;

cellnumGrowth = zeros(runNum,Nsp);
cellnumSel = zeros(runNum,Nsp);
TArr = cell(runNum,1);
NonLinFlagArr = zeros(runNum,Nsp);

lambdaSim={};
AgCapPopSim = {};
AgCapMeanSim = {};
% AgCapMeanEpitopesSim = {};
AgCapSTDSim = {};
EAgAbSim = {};
EAgAbExitSim = {};
AgCapExitSim = {};
BirthNumSim = {};
MutDistSim = {};

SeedingCellsSim = {};

OutEngsCellSim = {};


qSim = {};
q1Sim = {};
q2Sim = {};
PopulationTrack={};
ICArrtSim = {};
% Replace for 2224!
% flag_str = [NP1s,'_',NP2s,'_Nsp',num2str(Nsp),'_fAd',num2str(FitAdv),'_OnRCF',num2str(OnTargetRCF),'_OffRCF',num2str(OffTargetRCF)];
flag_str = [NP1s,'_',NP2s,'_Nsp',num2str(Nsp),'_fAd',num2str(FitAdv),'_OnRCF',num2str(OnTargetRCF),'_OffRCF',num2str(OffTargetRCF),'_id',num2str(id)];
runNum

aff = struct('BondsNum',BondsNum,'NBCR',NBCR,'SecBind',SecBind,'paramFlag',paramFlag,'EAbAg',EAbAg,'EAgMem',EAgMem,'qMFPTImmuneComp',qMFPTImmuneComp,'q1',q1,'q2',q2,...
    'a',a,'A',A,'NAg',NAg,'BindFlag',BindFlag,'F',F,'qMutP',qMutP,'tau',tau,'ClusterSize',ClusterSize,'CT',CT,'wcv',wcv,'rcv',rcv,'CorrVC',CorrVC,...
    'CorrVV',CorrVV,'MutVC',MutVC,'MutVV',MutVV,'ExtcFlag',ExtcFlag,'Epitope',1,'KonArr',KonArr,'AgThExp',AgThExp);

BCR_s = struct('l',lambda,'mu',mu,'D',Diff,'n',0,'dt',dt,'ParentNum',0,'Jackpotflag',Jackpotflag,'JPThreshold',JPThreshold,...
    'Nmax',Nmax,'lambdaBase',lambda0,'lambdaBaseExp',lambda0Exp,'aff',aff,'b',bOff,'muExpDecay',muExpDecay,'Seednum',0,'Origin',0)
%b stands for "boost" recieved per clone from interaction with HA
AgCap = [];
EAgAb = [];
AgCap_dist = [];
H = [];
for k = 1:runNum
    TMemoryGrowth = T1+T0;
    k
    OnTargetClones = OnTargetClonesNPSS; % NPHA
    
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    
%     ICArrt = {};
    ICArrt = cell(2,(T1+T2+T3+T4)/dt+4);
    
    cellNumDynComp = zeros((T1+T2+T3+T4)/dt+4,Nsp);
    cellNumDynGrowth = zeros((4*T0)/dt+4,Nsp);
    AgCapMean = zeros((T1+T2+T3+T4)/dt+4,1);
    AgCapSTD = zeros((T1+T2+T3+T4)/dt+4,1);
    lambdaPop = cell((T1+T2+T3+T4)/dt+4,Nsp);
    AgCapPop = cell((T1+T2+T3+T4)/dt+4,Nsp);   
    EAgAbPop = cell((T1+T2+T3+T4)/dt+4,Nsp);   
    AgICprop = [1 0 ; 1 0]; %%%%%%%%%remove
    global Virus;
    Virus= AgICprop(1,1);  
    global Virustype
      Virustype=1;
    ICPrepare(NAg,A,AgICprop,10*ICTypeNums,NBCR);

%     lambdaPop = cell(1,Nsp);
%     EAgAbPop = cell(1,Nsp);
    EAgAbExitPop = cell(2,Nsp);
    AgCapExitPop = cell(2,Nsp);
    
%     AgCapPop = cell(1,Nsp);
    BirthNum = cell(1,Nsp);
    MutDist = cell(4,Nsp);
    
    SeedingCells = cell(4,1);
    
    qPop = {};
    q1Pop = {};
    q2Pop = {};
    AgCapMean = [];
    AgCapSTD = [];
    AgCapMeanEpitopes = [];
    NspArr = [1:1:Nsp];
    exitFlag = 0;
    FlagFirst = 0;
    
    AgCapEpitopeCell = {};
    
    OutEngsRun = {};
    
    clear CellsArrSp;
    
    %%%%%%%% Remove B cells during the GCR %%%
    CellsArrSpExit = {};
    %%%%%%%%
    CellCount = 1;
    for n=1:Nsp
        clear bcArrClone;
        bcArrClone = Bcell.empty;
        for m=1:CellsPerCloneArr(n)
            
            % half of the cells have good affinity to Ag1 and the other half to
            % Ag2. They have bad affinity to the other type.
            BCR_s1 = BCR_s;
            
            BCR_s1.aff.EAbAg(1) = 0;
            
            %         BCR_s1.aff.EAbAg(2) = 5;
            
            if(find(OnTargetClones==n))
                BCR_s1.aff.EAbAg(2) = (-1+2*rand)*Multiplier +FitAdv;
	%	BCR_s1.aff.EAbAg(2) =-.25
                BCR_s1.aff.CorrVC =0;
                BCR_s1.aff.CorrVV =ones(length(CorrVV));
                BCR_s1.aff.MutVC =0;
                BCR_s1.aff.MutVV =0;
                BCR_s1.aff.EAbAg(3) = BCR_s1.aff.EAbAg(2);
                BCR_s1.aff.EAbAg(4) = BCR_s1.aff.EAbAg(2);%same epitope
            else
                 copulavar=(-1+2*(copularnd('t',CorrVV,1,1)))*Multiplier+FitAdv; %copula generates random numbers from the same distribution with given corr

                BCR_s1.aff.EAbAg(2) = copulavar(1);
                BCR_s1.aff.EAbAg(3) = copulavar(2);
                BCR_s1.aff.EAbAg(4) = copulavar(3);
	%	BCR_s1.aff.EAbAg(2) =-.25
                BCR_s1.aff.CorrVC =0;
                BCR_s1.aff.CorrVV =CorrVV;
                BCR_s1.aff.MutVC =0;
                BCR_s1.aff.MutVV =0;
                     end
            
 
            
            % Determine the epitope
            
            BCR_s1.aff.Epitope = n;
            BCR_s1.b = bOff;
            BCR_s1.aff.CT = COff*BCR_s1.aff.CT
            BCR_s1.Seednum=0;
            BCR_s1.Origin=0;
% 'Seednum',0,'Origin',0)
            
            bcArrClone(m,1) = Bcell(BCR_s1);
            
        end
        CellsArrSp{n} = bcArrClone;
    end

    for n=1:228
                    
        clear bcArrClone;
        bcArrClone = Bcell.empty;

        if(find(OnTargetClones==n))

            for m=1:memnumber
            
            % half of the cells have good affinity to Ag1 and the other half to
            % Ag2. They have bad affinity to the other type.
            BCR_s1 = BCR_s;
            
            BCR_s1.aff.EAbAg(1) = 0;
             BCR_s1.aff.EAbAg(2) = (-1 +2*rand)*Multiplier+FitAdv;
	%	BCR_s1.aff.EAbAg(2) =-.25
                BCR_s1.aff.CorrVC =0;
                BCR_s1.aff.CorrVV =ones(length(CorrVV));
                BCR_s1.aff.MutVC =0;
                BCR_s1.aff.MutVV =0;
                BCR_s1.aff.EAbAg(3) = BCR_s1.aff.EAbAg(2); %same epitope
            BCR_s1.aff.EAbAg(4) = BCR_s1.aff.EAbAg(2);
            % Determine the epitope
            
            BCR_s1.aff.Epitope = n;
            BCR_s1.b = bOff;
            BCR_s1.aff.CT = COff*BCR_s1.aff.CT;
            BCR_s1.Seednum=0;
            BCR_s1.Origin=0;
            
            bcArrClone(m,1) = Bcell(BCR_s1);
            bcArrClone(m,1).Mature=1;
            end
           CellsArrSp{n} = [CellsArrSp{n}; bcArrClone];
        else
            if CellsPerCloneArr(n)==0
            else
           for m=[] %change this back later.
                BCR_s1 = BCR_s;
            
            BCR_s1.aff.EAbAg(1) = 0;
%             copulavar=(-1+2*(copularnd('t',CorrVV,1,1)))*Multiplier+FitAdv; %copula generates random numbers from the same distribution with given corr
            copulavar=(-1+2*(copularnd('t',CorrVV,1,1)))*Multiplier+FitAdv; %copula generates random numbers from the same distribution with given corr
     
                BCR_s1.aff.EAbAg(2) = copulavar(1);
                BCR_s1.aff.EAbAg(3) = copulavar(2);
                BCR_s1.aff.EAbAg(4) = copulavar(3);
                BCR_s1.aff.CorrVC =0;
                BCR_s1.aff.CorrVV =CorrVV;
                BCR_s1.aff.MutVC =0;
                BCR_s1.aff.MutVV =0;
                          BCR_s1.aff.Epitope = n;
            BCR_s1.b = bOff;
            BCR_s1.aff.CT = COff*BCR_s1.aff.CT;
            BCR_s1.Seednum=0;
            BCR_s1.Origin=0;
            
            bcArrClone(m,1) = Bcell(BCR_s1);
            bcArrClone(m,1).Mature=1;
            end
  
           end
           CellsArrSp{n} = [CellsArrSp{n}; bcArrClone];
            end
        end
    

    MuFlag = 0;
    
%     cellNumDyn = zeros(1,Nsp);
    C = 1;
    
%     MemoryCellsArrSp = cell(1,Nsp);
    OutSideCellsArrSp = cell(1,Nsp);
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    %%% G1 %%%%%%%%
    %
    % Select cells to seed the initial GC (B1) Under NPHA - proxy for
    % Prime
%     ExitProb_InitalSeed=1;
%     lambdaMeanAll = 0;
%     CellsArrSpExit = AccumulateExitBCs(CellsArrSpExit,CellsArrSp,Nsp,ExitProb_InitalSeed,lambdaMeanAll);
%     res = SelectSeedAgCapTh(CellsArrSpExit,EAgAbExitPop,Nsp,paramFlag,1,1,CT/2);
%     %     res = SelectSeedAgCapThGermLine2(CellsArrSpExit,EAgAbExitPop,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,CT/2,...
%     %         FitAdv,OnTargetClonesNPSS,AllClonesNPSS,BCR_s,bOff,COff,CellsPerCloneArr);
%     CellsArrSp = res{1};
%     %     res = SelectSeedCloneNspNum(CellsArrSp,NumberOfCells4Seed);
%     res = SelectSeedCloneNspNumMemory(CellsArrSp,NumberOfCells4Seed);
%     
%     CellsArrSp = res{1};
%     NspArr = res{2};
%     OutSideCellsArrSp = res{3};
%     CellsArrSpExit = {};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CellsArrSpExit = {};

    AgDepClones = AgDepClonesList(NP);
%     AgThresholdGCSeed = ImmunogenicityAgTh(NP);
    AgThresholdGCSeed = ImmunClass.ImmunogenAgTh(NP);
AgThresholdGCSeed=100;
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSp,Nsp);

    displayCellNum(OutSideCellsArrSp(OnTargetClonesNPSS),'on-target stem prec');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPSS),'off-target stem prec');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPHA),'off-target head prec');
 
%     res =  SelectSeedAgCapThMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,CT/2,AgDepClones);
    res = SelectSeedAgCapThesholdMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,AgThresholdGCSeed,AgDepClones);

    CellsArrSp = res{1};
    OutSideCellsArrSp = res{3};
    res = SelectSeedCloneNspNumMemory(CellsArrSp,NumberOfCells4Seed);
    CellsArrSp = res{1};
    NspArr = res{2};
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{3},Nsp);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Which cells seed %%%
    
    SeedingCells{1,1} = SeedingCellsCheck(CellsArrSp,NspArr);
    
    OutEngsSeedPrime = CalcEngMatureNaive(CellsArrSp);
    OutEngsRun.OutEngsSeedPrime = OutEngsSeedPrime;
    displayCellNum(CellsArrSp(OnTargetClonesNPSS),'on-target stem Prime seed');
    displayCellNum(CellsArrSp(OffTargetClonesNPSS),'off-target stem Prime seed');
    displayCellNum(CellsArrSp(OffTargetClonesNPHA),'off-target head Prime seed');
    
    %%%%%%%%
    
    
    displayCellNum(OutSideCellsArrSp,'Outside cells Prime');
    
    displayCellNum(CellsArrSp,'Seed Prime');  
    
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%
    
    Cells2Seed = zeros(1,Nsp);
    for i=1:length(NspArr)
        n = NspArr(i);
        Cells2Seed(k,n) = length(CellsArrSp{n});
    end
        ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
 
    ICPrepare(NAg,A,AgICprop,ICTypeNums,NBCR);
    
%     res = GrowthPhase(CellsArrSp,NspArr,cellNumDyn,C,dt,T0,Nsp);
    C=1;
    res = GrowthPhase(CellsArrSp,NspArr,cellNumDynGrowth,C,dt,T0,Nsp);
    CellsArrSp = res{1};
    NspArr = res{2};
%     cellNumDyn = res{3};
    cellNumDynGrowth = res{3};
    C = res{4};
    
    
    C = 1;
    Pop = find(cellNumDynGrowth(C,:));
    SeedArr1 = SeedingCells{1}(:,1)';
    NspArr;
    if(length(setdiff(Pop,SeedArr1))>0)
        find(Cells2Seed)
    end
    
    
    displayCellNum(CellsArrSp,'Growth Prime');
    
%     PopulationTrack{k,1}=single(cellNumDyn);
    PopulationTrack{k,1}=single(cellNumDynGrowth);
    for i=1:length(NspArr)
        n = NspArr(i);
        cellnumGrowth(k,n) = length(CellsArrSp{n});
    end
    
    %%%%%%%%%%%
    %%%%%%%%%%
    
%     if(length(NspArr)==0)
%         continue;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Prime %%%%%%%
    
    MuFlag = 1;

    qProlCells_t = {};
    Countq = 1;
    qProlCells = {};
    

    CellseedComp = zeros(1,Nsp);
    for i=1:length(Nsp)
        CellseedComp(n) = length(CellsArrSp{n});
    end
    find(CellseedComp);
    
    C = 1;

%     Res = CompetGC(dt,T1,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDyn,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    Res = CompetGC(dt,T1,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    ;
    C=Res{1};
    CellsArrSpExit=Res{2};
    CellsArrSp=Res{3};
    Nsp=Res{4};
    ExitProb=Res{5};
%     cellNumDyn=Res{6};
    cellNumDynComp=Res{6};
    AgCapMean=Res{7};
    AgCapSTD=Res{8};
    lambdaPop=Res{9};
    AgCapPop=Res{10};
    ICArrt=Res{11};
    NspArr=Res{12};
    
%        Evar{k}=Res{15};
%     EMean{k}=Res{13};
%     EMed{k}=Res{14};
%     keyboard;
for n=1:228
  for j= 1:length(CellsArrSpExit{1,n})
      CellsArrSpExit{1,n}(j,1).BCR_s.Origin=1;
  end
end
    OutEngsExitPrime = CalcEngMatureNaive(CellsArrSpExit);
    OutEngsRun.OutEngsExitPrime = OutEngsExitPrime;
    displayCellNum(CellsArrSpExit(OnTargetClonesNPSS),'on-target stem Prime exit');
    displayCellNum(CellsArrSpExit(OffTargetClonesNPSS),'off-target stem Prime exit');
    displayCellNum(CellsArrSpExit(OffTargetClonesNPHA),'off-target head Prime Exit');
    
    CellAfterComp = zeros(1,Nsp);
    for i=1:length(Nsp)
        CellAfterComp(n) = length(CellsArrSp{n});
    end
    find(CellAfterComp);
    
    
    C = 1;
    Pop1 = find(cellNumDynComp(C,:));
    
% % % %     OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSpExit,Nsp);
   

%     displayCellNum(OutSideCellsArrSp,'Total Prime');
%     displayCellNum(CellsArrSpExit,'Exit Prime');
%     disp([num2str(NN),' Mature B cells Prime']);
%     displayCellNum(CellsArrSp,'Comp Prime');
%     
%     displayCellNum(OutSideCellsArrSp(OnTargetClonesNPSS),'on-target stem Prime');
%     displayCellNum(OutSideCellsArrSp(OffTargetClonesNPSS),'off-target stem Prime');
%     displayCellNum(OutSideCellsArrSp(OffTargetClonesNPHA),'off-target head Prime');
    
% %     OutEngsAfterPrime = CalcEngMatureNaive(OutSideCellsArrSp);
% %     OutEngsRun.OutEngsAfterPrime = OutEngsAfterPrime;
    
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Expand memory B cells Prime %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     We simulate the Ag dependent proliferation of memory B cell as a
%     "fake" GC cycle of TMemoryGrowth days

    %%% Initiate ICs %%%
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
%     AgICprop = [1 0 ; 1 0]
    ExitProbOUtSideGC = 0; 

    ICPrepare(NAg_T2,A,AgICprop,ICMemExpMultiplier*ICTypeNums,NBCR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Amongst the cells that were left out and did not enter the GC we
    %%%% select the memory cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AgThresholdGCSeed
%     res = SplitMemoryNaive(OutSideCellsArrSp,AgDepClones,Nsp);
    res = SplitMemoryNaive2(OutSideCellsArrSp,AgDepClones,Nsp);
%     res = SplitMemoryNaiveAgDep(OutSideCellsArrSp,AgDepClones,Nsp,AgThresholdGCSeed);
    CellsArrMemorySp = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};

    % Decay of Ag independent memory B cells
    CellArrMemAgInDep = res{4}; % Ag independent memory B cell will decay
    NspArrMemAgInDep = res{5};
    displayCellNum(CellArrMemAgInDep,'Number of Ag independent memory B cells Prime');
    res = ExpansionDecay(CellArrMemAgInDep,NspArrMemAgInDep,cellNumDynGrowth,C,dt,TExpansionDecayAgID,Nsp); % Decay of the Ag ID
    displayCellNum(res{1},'Number of Ag independent memory B cells Prime after decay');
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{1},Nsp);

    res=SplitMemoryWeak(CellsArrMemorySp,AgDepClones,Nsp,Agcutoff);
             CellsArrMemorySp = res{1};
%     NspArrMemory = res{2};
    CellsArrMemorySpNonVia = res{3}; 
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{3},Nsp);
    
    OutEngsExpPreGrowthPrime = CalcEngMatureNaive(CellsArrMemorySp);
    OutEngsRun.OutEngsExpPreGrowthPrime = OutEngsExpPreGrowthPrime;
    displayCellNum(CellsArrMemorySp,'Number of Ag specific memory B cells Prime');
    % Initial memory B cell proliferation
    CellsArrSpMemoryExit = {};


    Res = GrownMemBcells(dt,TMemoryGrowth,NspArrMemory,C,CellsArrSpMemoryExit,CellsArrMemorySp,Nsp,ExitProbOUtSideGC,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt,AgThExp);
   %     Res = GrownMemBcells(dt,TMemoryGrowth,NspArrMemory,C,CellsArrSpMemoryExit,CellsArrMemorySp,Nsp,ExitProbOUtSideGC,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt,AgThExp);
    CellsArrMemorySpProf=Res{3};
    for n=1:228
        for j= 1:length(CellsArrMemorySpProf{1,n})
            CellsArrMemorySpProf{1,n}(j,1).BCR_s.Origin=2;
        end
    end
    OutEngsExpPostGrowthPrime = CalcEngMatureNaive(CellsArrMemorySpProf);
    OutEngsRun.OutEngsExpPostGrowthPrime = OutEngsExpPostGrowthPrime;
    displayCellNum(CellsArrMemorySpProf,'Number of Ag specific memory B cells Prime following proliferation');
    NspArrMemoryProf = Res{12};
    
    displayCellNum(CellsArrMemorySpProf(OnTargetClonesNPSS),'on-target stem prime Expansion');
    displayCellNum(CellsArrMemorySpProf(OffTargetClonesNPSS),'off-target stem prime Expansion');
    displayCellNum(CellsArrMemorySpProf(OffTargetClonesNPHA),'off-target head prime Expansion');
    
    % Merge the B cells created in the GC with those outside
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSpExit,Nsp); 
    
    displayCellNum(OutSideCellsArrSp,'Prime: Naive and mem created in the GC. Before expansion');
    
    % Merge the proliferated B cell population with the population outside
  %   OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    OutSideCellsArrSpAllCells = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    
    
NN = 0;
    MM = 0;
    if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
        for n=1:Nsp
            if(~isempty(OutSideCellsArrSp{n}))
                MM = MM + length(OutSideCellsArrSp{n});
                idx = find([OutSideCellsArrSp{n}.Mature]);
                if(length(idx))
                    NN = NN + length(idx);
                    EAgAbExitPop{1,n} = reshape(GetAffField(OutSideCellsArrSp{n}(idx),'EAbAg'),NAnt+1,length(OutSideCellsArrSp{n}(idx)))';
                    AgCapExitPop{1,n} = [OutSideCellsArrSp{n}(idx).AgCap];
%                     MutDist{1,n} = [OutSideCellsArrSp{n}.ParentNum];
                    BCRsCells = [OutSideCellsArrSp{n}.BCR_s];
                    MutDist{1,n} = [BCRsCells.ParentNum];
                end
            end
        end
    end
    
    displayCellNum(OutSideCellsArrSpAllCells,'Total Prime. After expansion');
    displayCellNum(CellsArrSpExit,'Exit Prime');
    disp([num2str(NN),' Mature B cells Prime']);
    displayCellNum(CellsArrSp,'Comp Prime');
 
%finding ALL Memory.
res = SplitMemoryNaive2(OutSideCellsArrSpAllCells,AgDepClones,Nsp);
%     res = SplitMemoryNaiveAgDep(OutSideCellsArrSp,AgDepClones,Nsp,AgThresholdGCSeed);
    CellsArrMemorySpProf = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};

    %%%%% The expanded populaiton decayed %%%%%%%   
%     TExpansionDecay = TMemoryGrowth;
    res = ExpansionDecay(CellsArrMemorySpProf,NspArrMemory ,cellNumDynGrowth,C,dt,TExpansionDecay,Nsp);
    CellsArrMemorySpProfLeft = res{1};
    displayCellNum(CellsArrMemorySpProfLeft,'Expansion decay memory B cells left Prime');
    
    %OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProfLeft,Nsp);
OutSideCellsArrSp =CellsArrMemorySpProfLeft;

    displayCellNum(OutSideCellsArrSp,'Total Prime. After expansion decay');
    
    
    displayCellNum(OutSideCellsArrSp(OnTargetClonesNPSS),'on-target stem Prime after Expansion decay');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPSS),'off-target stem Prime after Expansion decay');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPHA),'off-target head Prime after Expansion decay');
    
    
    %%%%%%%%%%%%%%% B1 %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     CellsArrSpExit = {};
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    AgICprop = [0 1 ; 0 1];
    Virus= AgICprop(1,1);
    ICPrepare(NAg_T2,A,AgICprop,ICTypeNums,NBCR);
%     NumberOfCells4Seed = 100    
    NP = NP2s;
%     OnTargetClones = [11 13 14 152 167 169 170 171 172 173]; %NPSS
    OnTargetClones = ClonesList(NP);
    AgDepClones = AgDepClonesList(NP);
%     AgThresholdGCSeed = ImmunogenicityAgTh(NP);
%     AgThresholdGCSeed = ImmunClass.ImmunogenAgTh(NP);
    
%     KonArr = KonEpitopes(NP);
    KonArr = KonEpitopes(NP,q01Arm,q02Arm);
    
    aff = struct('BondsNum',BondsNum,'NBCR',NBCR,'SecBind',SecBind,'paramFlag',paramFlag,'EAbAg',EAbAg,'EAgMem',EAgMem,'qMFPTImmuneComp',qMFPTImmuneComp,'q1',q1,'q2',q2,...
        'a',a,'A',A,'NAg',NAg,'BindFlag',BindFlag,'F',F,'qMutP',qMutP,'tau',tau,'ClusterSize',ClusterSize,'CT',CT,'wcv',wcv,'rcv',rcv,'CorrVC',CorrVC,...
        'CorrVV',CorrVV,'MutVC',MutVC,'MutVV',MutVV,'ExtcFlag',ExtcFlag,'Epitope',1,'KonArr',KonArr,'AgThExp',AgThExp);
    
    BCR_s = struct('l',lambda,'mu',mu,'D',Diff,'n',0,'dt',dt,'ParentNum',0,'Jackpotflag',Jackpotflag,'JPThreshold',JPThreshold,...
        'Nmax',Nmax,'lambdaBase',lambda0,'lambdaBaseExp',lambda0Exp,'aff',aff,'b',bOff,'muExpDecay',muExpDecay);  
      
    for n=1:Nsp
        if(~isempty(OutSideCellsArrSp{n}))
            for j=1:length(OutSideCellsArrSp{n})
                OutSideCellsArrSp{n}(j).BCR_s.aff.KonArr = KonArr;
                OutSideCellsArrSp{n}(j).AgCap = OutSideCellsArrSp{n}(j).AffinityClass.CAgCaptureFunc();
            end
        end
    end
    
    OutEngsB1Init = CalcEngMatureNaive(OutSideCellsArrSp);
    OutEngsRun.OutEngsB1Init = OutEngsB1Init;
    
    % Update kOn for all B cells, and Ag capture
    
    
    %%% Select cells to seed %%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     res = SelectSeedAgCapThMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,CT/2,AgDepClones);
    res = SelectSeedAgCapThesholdMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,AgThresholdGCSeed,AgDepClones);
    CellsArrSp = res{1};
    OutSideCellsArrSp = res{3};
    res = SelectSeedCloneNspNumMemory(CellsArrSp,NumberOfCells4Seed);
    CellsArrSp = res{1}; % cells selected to enter GC
    NspArr = res{2};
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{3},Nsp); % the cells left outside of the GC
    
    %%% Which cells seed %%%
    
    SeedingCells{2,1} = SeedingCellsCheck(CellsArrSp,NspArr);
    
    %%%%%%%%
    
    displayCellNum(CellsArrSp(OnTargetClonesNPSS),'on-target stem B1 seed');
    displayCellNum(CellsArrSp(OffTargetClonesNPSS),'off-target stem B1 seed');
    displayCellNum(CellsArrSp(OffTargetClonesNPHA),'off-target head B1 seed');
        CellsArrSp
        seedcount=1;
        for n=1:228
  for j= 1:length(CellsArrSp{1,n})
      CellsArrSp{1,n}(j,1).BCR_s.Seednum=seedcount;
      seedcount=seedcount+1;
  end
end
    OutEngsSeedB1 = CalcEngMatureNaive(CellsArrSp);
    OutEngsRun.OutEngsSeedB1 = OutEngsSeedB1;
    displayCellNum(CellsArrSp,'Seed B1');
    displayCellNum(OutSideCellsArrSp,'Outside cells B1');

    %%%%%% B1 %%%%%%%%%

    C = T0/dt +2;
    res = GrowthPhase(CellsArrSp,NspArr,cellNumDynGrowth,C,dt,T0,Nsp);
    CellsArrSp = res{1};
    NspArr = res{2};
%     cellNumDyn = res{3};
    cellNumDynGrowth = res{3};
    C = res{4};
    
    displayCellNum(CellsArrSp,'Growth B1');
    
    %%%%%%%%
    CellsArrSpExit = {};
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    
    ICPrepare(NAg_T2,A,AgICprop,ICTypeNums,NBCR);
    
    %%%%%%%%
    C = T1/dt +2;
    
    %%%%%%%%

    Res = CompetGC(dt,T2,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    C=Res{1};
    CellsArrSpExit=Res{2};
    CellsArrSp=Res{3};
    Nsp=Res{4};
    ExitProb=Res{5};
%     cellNumDyn=Res{6};
    cellNumDynComp=Res{6};
    AgCapMean=Res{7};
    AgCapSTD=Res{8};
    lambdaPop=Res{9};
    AgCapPop=Res{10};
    ICArrt=Res{11};
    NspArr=Res{12};
  
    OutEngsExitB1 = CalcEngMatureNaive(CellsArrSpExit);
    OutEngsRun.OutEngsExitB1 = OutEngsExitB1;
    displayCellNum(CellsArrSpExit(OnTargetClonesNPSS),'on-target stem B1 exit');
    displayCellNum(CellsArrSpExit(OffTargetClonesNPSS),'off-target stem B1 exit');
    displayCellNum(CellsArrSpExit(OffTargetClonesNPHA),'off-target head B1 Exit');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Expand memory B cells B1 %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     We simulate the Ag dependent proliferation of memory B cell as a
%     "fake" GC cycle of TMemoryGrowth days

    %%% Initiate ICs %%%
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
%     AgICprop = [1 0 ; 1 0]
    ExitProbOUtSideGC = 0; 

    ICPrepare(NAg_T2,A,AgICprop,ICMemExpMultiplier*ICTypeNums,NBCR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Amongst the cells that were left out and did not enter the GC we
    %%%% select the memory cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AgThresholdGCSeed
%     res = SplitMemoryNaive(OutSideCellsArrSp,AgDepClones,Nsp);
    res = SplitMemoryNaive2(OutSideCellsArrSp,AgDepClones,Nsp);
%     res = SplitMemoryNaiveAgDep(OutSideCellsArrSp,AgDepClones,Nsp,AgThresholdGCSeed);
    CellsArrMemorySp = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};
    
    % Decay of Ag independent memory B cells
    CellArrMemAgInDep = res{4}; % Ag independent memory B cell will decay
    NspArrMemAgInDep = res{5};
    displayCellNum(CellArrMemAgInDep,'Number of Ag independent memory B cells B1');
    res = ExpansionDecay(CellArrMemAgInDep,NspArrMemAgInDep,cellNumDynGrowth,C,dt,TExpansionDecayAgID,Nsp); % Decay of the Ag ID
    displayCellNum(res{1},'Number of Ag independent memory B cells B1 after decay');
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{1},Nsp);
    
     res=SplitMemoryWeak(CellsArrMemorySp,AgDepClones,Nsp,Agcutoff);
             CellsArrMemorySp = res{1};
%     NspArrMemory = res{2};
    CellsArrMemorySpNonVia = res{3}; 
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{3},Nsp);

         for n=1:228
        for j= 1:length(CellsArrMemorySpProf{1,n})
            CellsArrMemorySpProf{1,n}(j,1).BCR_s.Origin=2;
        end
    end


    OutEngsExpPreGrowthB1 = CalcEngMatureNaive(CellsArrMemorySp);
    OutEngsRun.OutEngsExpPreGrowthB1 = OutEngsExpPreGrowthB1;
    displayCellNum(CellsArrMemorySp,'Number of Ag specific memory B cells B1');
    % Initial memory B cell proliferation
    CellsArrSpMemoryExit = {};

    Res = GrownMemBcells(dt,TMemoryGrowth,NspArrMemory,C,CellsArrSpMemoryExit,CellsArrMemorySp,Nsp,ExitProbOUtSideGC,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt,AgThExp);
    CellsArrMemorySpProf=Res{3};
    TMemoryGrowth=0; %only if you don't want an  SPF
    OutEngsExpPostGrowthB1 = CalcEngMatureNaive(CellsArrMemorySpProf);
    OutEngsRun.OutEngsExpPostGrowthB1 = OutEngsExpPostGrowthB1;
    displayCellNum(CellsArrMemorySpProf,'Number of Ag specific memory B cells B1 following proliferation');
    NspArrMemoryProf = Res{12};
    
    displayCellNum(CellsArrMemorySpProf(OnTargetClonesNPSS),'on-target stem B1 Expansion');
    displayCellNum(CellsArrMemorySpProf(OffTargetClonesNPSS),'off-target stem B1 Expansion');
    displayCellNum(CellsArrMemorySpProf(OffTargetClonesNPHA),'off-target head B1 Expansion');
    
    % Merge the B cells created in the GC with those outside
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSpExit,Nsp); 
    
    displayCellNum(OutSideCellsArrSp,'B1: Naive and mem created in the GC. Before expansion');
    
    % Merge the proliferated B cell population with the population outside
%     OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    OutSideCellsArrSpAllCells = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    
    NN = 0;
    if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
        for n=1:Nsp
            if(~isempty(OutSideCellsArrSpAllCells{n}))
                idx = find([OutSideCellsArrSpAllCells{n}.Mature]);
                if(length(idx))
                    NN = NN + length(idx);
                    EAgAbExitPop{2,n} = reshape(GetAffField(OutSideCellsArrSpAllCells{n}(idx),'EAbAg'),NAnt+1,length(OutSideCellsArrSpAllCells{n}(idx)))';
                    AgCapExitPop{2,n} = [OutSideCellsArrSpAllCells{n}(idx).AgCap];
%                     MutDist{2,n} = [OutSideCellsArrSpAllCells{n}.ParentNum];
                    BCRsCells = [OutSideCellsArrSpAllCells{n}.BCR_s];
                    MutDist{2,n} = [BCRsCells.ParentNum];
                end
            end
        end
    end
    
    %%%%%%%%
%     if(length(CellsArrSpExit))
%         if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
%             for n=1:Nsp
%                 if(~isempty(CellsArrSpExit{n}))
%                     EAgAbExitPop{2,n} = reshape(GetAffField(CellsArrSpExit{n},'EAbAg'),3,length(CellsArrSpExit{n}))';
%                     AgCapExitPop{2,n} = [CellsArrSpExit{n}.AgCap];
%                 end
%             end
%         end
%     end

    displayCellNum(OutSideCellsArrSpAllCells,'Total B1. After expansion');
    displayCellNum(CellsArrSpExit,'Exit B1');
    disp([num2str(NN),' Mature B cells B1']);
    displayCellNum(CellsArrSp,'Comp B1');

 res = SplitMemoryNaive2(OutSideCellsArrSpAllCells,AgDepClones,Nsp);
%     res = SplitMemoryNaiveAgDep(OutSideCellsArrSp,AgDepClones,Nsp,AgThresholdGCSeed);
    CellsArrMemorySpProf = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};

    %%%%% The expanded populaiton decayed %%%%%%%   
%     TExpansionDecay = TMemoryGrowth;
    res = ExpansionDecay(CellsArrMemorySpProf,NspArrMemory ,cellNumDynGrowth,C,dt,TExpansionDecay,Nsp);
    CellsArrMemorySpProfLeft = res{1};
    displayCellNum(CellsArrMemorySpProfLeft,'Expansion decay memory B cells left B1');
    
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProfLeft,Nsp);
    displayCellNum(OutSideCellsArrSp,'Total B1. After expansion decay');
    
    
    displayCellNum(OutSideCellsArrSp(OnTargetClonesNPSS),'on-target stem B1 after Expansion decay');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPSS),'off-target stem B1 after Expansion decay');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPHA),'off-target head B1 after Expansion decay');
    
    
    %%%%%%%%%%% Update Ag capture before seedning%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
%     AgICprop = [1 0 ; 1 0]
%     Virus= AgICprop(1,1);
    ICPrepare(NAg_T2,A,AgICprop,ICTypeNums,NBCR);
    NP = NP3s;
    OnTargetClones = ClonesList(NP);
    AgDepClones = AgDepClonesList(NP);
%     AgThresholdGCSeed = ImmunogenicityAgTh(NP);
%     AgThresholdGCSeed = ImmunClass.ImmunogenAgTh(NP);
    
%     KonArr = KonEpitopes(NP);
    KonArr = KonEpitopes(NP,q01Arm,q02Arm);
    
    aff = struct('BondsNum',BondsNum,'NBCR',NBCR,'SecBind',SecBind,'paramFlag',paramFlag,'EAbAg',EAbAg,'EAgMem',EAgMem,'qMFPTImmuneComp',qMFPTImmuneComp,'q1',q1,'q2',q2,...
        'a',a,'A',A,'NAg',NAg,'BindFlag',BindFlag,'F',F,'qMutP',qMutP,'tau',tau,'ClusterSize',ClusterSize,'CT',CT,'wcv',wcv,'rcv',rcv,'CorrVC',CorrVC,...
        'CorrVV',CorrVV,'MutVC',MutVC,'MutVV',MutVV,'ExtcFlag',ExtcFlag,'Epitope',1,'KonArr',KonArr,'AgThExp',AgThExp);
    
    BCR_s = struct('l',lambda,'mu',mu,'D',Diff,'n',0,'dt',dt,'ParentNum',0,'Jackpotflag',Jackpotflag,'JPThreshold',JPThreshold,...
        'Nmax',Nmax,'lambdaBase',lambda0,'lambdaBaseExp',lambda0Exp,'aff',aff,'b',bOff,'muExpDecay',muExpDecay);  
  
      Virustype=3;

    for n=1:Nsp
        if(~isempty(OutSideCellsArrSp{n}))
            for j=1:length(OutSideCellsArrSp{n})
                OutSideCellsArrSp{n}(j).BCR_s.aff.KonArr = KonArr;
                OutSideCellsArrSp{n}(j).AgCap = OutSideCellsArrSp{n}(j).AffinityClass.CAgCaptureFunc();
            end
        end
    end
    
    OutEngsB2Init = CalcEngMatureNaive(OutSideCellsArrSp);
    OutEngsRun.OutEngsB2Init = OutEngsB2Init;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% Second Boost %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%    
    
%     res = SelectSeedAgCapThMemory(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,CT/2);  
%     res = SelectSeedAgCapThMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,CT/2,AllClonesNPSS);
%     CellsArrSp = res{1};
%     OutSideCellsArrSp = res{3};
%     
% %     res = SelectSeedCloneNspNum(CellsArrSp,NumberOfCells4Seed);
% %     CellsArrSp = res{1};
% %     NspArr = res{2};
%     
%     res = SelectSeedCloneNspNumMemory(CellsArrSp,NumberOfCells4Seed);
%     CellsArrSp = res{1};
%     NspArr = res{2};

%     res = SelectSeedAgCapThMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,CT/2,AgDepClones);
    res = SelectSeedAgCapThesholdMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,AgThresholdGCSeed,AgDepClones);
    CellsArrSp = res{1};
    OutSideCellsArrSp = res{3};
    
    res = SelectSeedCloneNspNumMemory(CellsArrSp,NumberOfCells4Seed);
    CellsArrSp = res{1};
    NspArr = res{2};
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{3},Nsp);
    
    %%% Which cells seed %%%
    
    OutEngsSeedB2 = CalcEngMatureNaive(CellsArrSp);
    OutEngsRun.OutEngsSeedB2 = OutEngsSeedB2;
    SeedingCells{3,1} = SeedingCellsCheck(CellsArrSp,NspArr);
    
    %%%%%%%%

    displayCellNum(CellsArrSp,'Seed B2');
    
    displayCellNum(CellsArrSp(OnTargetClonesNPSS),'on-target stem B2 seed');
    displayCellNum(CellsArrSp(OffTargetClonesNPSS),'off-target stem B2 seed');
    displayCellNum(CellsArrSp(OffTargetClonesNPHA),'off-target head B2 seed');
    
    displayCellNum(OutSideCellsArrSp,'Outside cells B2');
    %%%%%% Growth Boost2 %%%%%%%%%
    
%     res = GrowthPhase(CellsArrSp,NspArr,cellNumDyn,C,dt,T0,Nsp);

    C = 2*T0/dt+3;
    res = GrowthPhase(CellsArrSp,NspArr,cellNumDynGrowth,C,dt,T0,Nsp);
    CellsArrSp = res{1};
    NspArr = res{2};
%     cellNumDyn = res{3};
    cellNumDynGrowth = res{3};
    C = res{4};
    
    displayCellNum(CellsArrSp,'Growth B2');
    
    %%%%%%%%
    CellsArrSpExit = {};
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    
    ICPrepare(NAg_T3,A,AgICprop,ICTypeNums,NBCR);
    
    %%%%%%%%
    
    C = T1/dt+T2/dt+3;
    
    %%%%%%%%
        
%     Res = CompetGC(dt,T3,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDyn,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    Res = CompetGC(dt,T3,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    C=Res{1};
    CellsArrSpExit=Res{2};
    CellsArrSp=Res{3};
    Nsp=Res{4};
    ExitProb=Res{5};
%     cellNumDyn=Res{6};
    cellNumDynComp=Res{6};
    AgCapMean=Res{7};
    AgCapSTD=Res{8};
    lambdaPop=Res{9};
    AgCapPop=Res{10};
    ICArrt=Res{11};
    NspArr=Res{12};
    
    OutEngsExitB2 = CalcEngMatureNaive(CellsArrSpExit);
    OutEngsRun.OutEngsExitB2 = OutEngsExitB2;
    displayCellNum(CellsArrSpExit(OnTargetClonesNPSS),'on-target stem B2 exit');
    displayCellNum(CellsArrSpExit(OffTargetClonesNPSS),'off-target stem B2 exit');
    displayCellNum(CellsArrSpExit(OffTargetClonesNPHA),'off-target head B2 Exit');
    
    %%%%%%%%
    %%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Expand memory B cells B2 %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     We simulate the Ag dependent proliferation of memory B cell as a
%     "fake" GC cycle of TMemoryGrowth days

    %%% Initiate ICs %%%
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    
    CellsArrSpMemoryExit = {};
    ExitProbOUtSideGC = 0; 
    ICPrepare(NAg_T2,A,AgICprop,ICMemExpMultiplier*ICTypeNums,NBCR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Amongst the cells that were left out and did not enter the GC we
    %%%% select the memory cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     res = SplitMemoryNaive(OutSideCellsArrSp,AgDepClones,Nsp);
    res = SplitMemoryNaive2(OutSideCellsArrSp,AgDepClones,Nsp);
    CellsArrMemorySp = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};
    
        % Decay of Ag independent memory B cells
    CellArrMemAgInDep = res{4}; % Ag independent memory B cell will decay
    NspArrMemAgInDep = res{5};
    displayCellNum(CellArrMemAgInDep,'Number of Ag independent memory B cells B2');
    res = ExpansionDecay(CellArrMemAgInDep,NspArrMemAgInDep,cellNumDynGrowth,C,dt,TExpansionDecayAgID,Nsp);
    displayCellNum(res{1},'Number of Ag independent memory B cells B2 after decay');
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{1},Nsp);
    
    OutEngsExpPreGrowthB2 = CalcEngMatureNaive(CellsArrMemorySp);
    OutEngsRun.OutEngsExpPreGrowthB2 = OutEngsExpPreGrowthB2;
    
    displayCellNum(CellsArrMemorySp,'Number of Ag specific memory B cells B2');
    % Initial memory B cell proliferation

    Res = GrownMemBcells(dt,TMemoryGrowth,NspArrMemory,C,CellsArrSpMemoryExit,CellsArrMemorySp,Nsp,ExitProbOUtSideGC,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt,AgThExp);
    CellsArrMemorySpProf=Res{3};
    displayCellNum(CellsArrMemorySpProf,'Number of Ag specific memory B cells B2 following proliferation');
    NspArrMemoryProf = Res{12};
    
    OutEngsExpPostGrowthB2 = CalcEngMatureNaive(CellsArrMemorySpProf);
    OutEngsRun.OutEngsExpPostGrowthB2 = OutEngsExpPostGrowthB2;
    
    displayCellNum(CellsArrMemorySpProf(OnTargetClonesNPSS),'on-target stem B2 Expansion');
    displayCellNum(CellsArrMemorySpProf(OffTargetClonesNPSS),'off-target stem B2 Expansion');
    displayCellNum(CellsArrMemorySpProf(OffTargetClonesNPHA),'off-target head B2 Expansion');
    
    % Merge the B cells created in the GC with those outside
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSpExit,Nsp); 
           
    displayCellNum(OutSideCellsArrSp,'B2: Naive and mem created in the GC. Before expansion');
    
%     Merge the proliferated B cell population with the population outside
%     OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    OutSideCellsArrSpAllCells = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    %%%%%%%%
    
    NN = 0;
    if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
        for n=1:Nsp
            if(~isempty(OutSideCellsArrSpAllCells{n}))
                idx = find([OutSideCellsArrSpAllCells{n}.Mature]);
                if(length(idx))
                    NN = NN + length(idx);
                    EAgAbExitPop{3,n} = reshape(GetAffField(OutSideCellsArrSpAllCells{n}(idx),'EAbAg'),NAnt+1,length(OutSideCellsArrSpAllCells{n}(idx)))';
                    AgCapExitPop{3,n} = [OutSideCellsArrSpAllCells{n}(idx).AgCap];
%                     MutDist{3,n} = [OutSideCellsArrSpAllCells{n}.ParentNum];
                    BCRsCells = [OutSideCellsArrSpAllCells{n}.BCR_s];
                    MutDist{3,n} = [BCRsCells.ParentNum];
                end
            end
        end
    end
    

    displayCellNum(OutSideCellsArrSpAllCells,'Total B2. After expansion');
    displayCellNum(CellsArrSpExit,'Exit B2');
    disp([num2str(NN),' Mature B cells B2']);
    displayCellNum(CellsArrSp,'Comp B2');


 res = SplitMemoryNaive2(OutSideCellsArrSpAllCells,AgDepClones,Nsp);
%     res = SplitMemoryNaiveAgDep(OutSideCellsArrSp,AgDepClones,Nsp,AgThresholdGCSeed);
    CellsArrMemorySpProf = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};


    %%%%% The expanded populaiton decayed %%%%%%%   
%     TExpansionDecay = TMemoryGrowth;
    res = ExpansionDecay(CellsArrMemorySpProf,NspArrMemoryProf,cellNumDynGrowth,C,dt,TExpansionDecay,Nsp);
    CellsArrMemorySpProfLeft = res{1};
    displayCellNum(CellsArrMemorySpProfLeft,'Expansion decay memory B cells left B2');
    
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProfLeft,Nsp);
    displayCellNum(OutSideCellsArrSp,'Total B2. After expansion decay');

    displayCellNum(OutSideCellsArrSp(OnTargetClonesNPSS),'on-target stem B2 after Expansion decay');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPSS),'off-target stem B2 after Expansion decay');
    displayCellNum(OutSideCellsArrSp(OffTargetClonesNPHA),'off-target head B2 after Expansion decay');
    
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Flu challenge! %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if(FCFlag)
    
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
%     AgICprop = [0 1 ; 0 1]
%     Virus= AgICprop(1,1);
    ICPrepare(NAg_T2,A,AgICprop,ICTypeNums,NBCR);
    
    NP = NP4s;
%     OnTargetClones = OnTargetClonesNPHA; %NPHA   
    OnTargetClones = ClonesList(NP);
%     AgThresholdGCSeed = ImmunogenicityAgTh(NP);
%     AgThresholdGCSeed = ImmunClass.ImmunogenAgTh(NP);
    AgDepClones = AgDepClonesList(NP);
    
%     KonArr = KonEpitopes(NP);
    KonArr = KonEpitopes(NP,q01Arm,q02Arm);
    
    aff = struct('BondsNum',BondsNum,'NBCR',NBCR,'SecBind',SecBind,'paramFlag',paramFlag,'EAbAg',EAbAg,'EAgMem',EAgMem,'qMFPTImmuneComp',qMFPTImmuneComp,'q1',q1,'q2',q2,...
        'a',a,'A',A,'NAg',NAg,'BindFlag',BindFlag,'F',F,'qMutP',qMutP,'tau',tau,'ClusterSize',ClusterSize,'CT',CT,'wcv',wcv,'rcv',rcv,'CorrVC',CorrVC,...
        'CorrVV',CorrVV,'MutVC',MutVC,'MutVV',MutVV,'ExtcFlag',ExtcFlag,'Epitope',1,'KonArr',KonArr,'AgThExp',AgThExp);
    
    BCR_s = struct('l',lambda,'mu',mu,'D',Diff,'n',0,'dt',dt,'ParentNum',0,'Jackpotflag',Jackpotflag,'JPThreshold',JPThreshold,...
        'Nmax',Nmax,'lambdaBase',lambda0,'lambdaBaseExp',lambda0Exp,'aff',aff,'b',bOff,'muExpDecay',muExpDecay);
    
    for n=1:Nsp
        if(~isempty(OutSideCellsArrSp{n}))
            for j=1:length(OutSideCellsArrSp{n})
                OutSideCellsArrSp{n}(j).BCR_s.aff.KonArr = KonArr;
                OutSideCellsArrSp{n}(j).AgCap = OutSideCellsArrSp{n}(j).AffinityClass.CAgCaptureFunc();
            end
        end
    end
    
        OutEngsInitChallenge = CalcEngMatureNaive(OutSideCellsArrSp);
    OutEngsRun.OutEngsInitChallenge = OutEngsInitChallenge;
    %%%%%%%%%%%%%%%%%%%%%%
    
    res = SelectSeedAgCapThesholdMemoryAgDep(OutSideCellsArrSp,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,AgThresholdGCSeed,AgDepClones);
    CellsArrSp = res{1};
    OutSideCellsArrSp = res{3};  
    res = SelectSeedCloneNspNumMemory(CellsArrSp,NumberOfCells4Seed);
    CellsArrSp = res{1};
    NspArr = res{2};
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{3},Nsp);
    
    %%% Which cells seed %%%
    
    SeedingCells{4,1} = SeedingCellsCheck(CellsArrSp,NspArr);
    
    %%%%%%%%

    displayCellNum(CellsArrSp,'Seed Flu challenge');
    displayCellNum(OutSideCellsArrSp,'Outside cells Flu challenge');
    %%%%%% Growth Boost2 %%%%%%%%%

    C = 3*T0/dt+4;
    res = GrowthPhase(CellsArrSp,NspArr,cellNumDynGrowth,C,dt,T0,Nsp);
    CellsArrSp = res{1};
    NspArr = res{2};
%     cellNumDyn = res{3};
    cellNumDynGrowth = res{3};
    C = res{4};
    
    displayCellNum(CellsArrSp,'Growth Flu challenge');
    
    
    %%%%%%%%
    CellsArrSpExit = {};
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    
    ICPrepare(NAg_T3,A,AgICprop,ICTypeNums,NBCR);
    
    %%%%%%%%
    
    C = T1/dt+T2/dt+T3/dt+4;
    
    %%%%%%%%
        
%     Res = CompetGC(dt,T4,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDyn,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    Res = CompetGC(dt,T4,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt);
    C=Res{1};
    CellsArrSpExit=Res{2};
    CellsArrSp=Res{3};
    Nsp=Res{4};
    ExitProb=Res{5};
%     cellNumDyn=Res{6};
    cellNumDynComp=Res{6};
    AgCapMean=Res{7};
    AgCapSTD=Res{8};
    lambdaPop=Res{9};
    AgCapPop=Res{10};
    ICArrt=Res{11};
    NspArr=Res{12};
    
    %%%%%%%%
    %%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Expand memory B cells Flu challenge %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     We simulate the Ag dependent proliferation of memory B cell as a
%     "fake" GC cycle of TMemoryGrowth days

    %%% Initiate ICs %%%
    ICArr = {};
    CI1 = cell(1e3);
    CI2 = cell(1e3);
    
    CellsArrSpMemoryExit = {};
    ExitProbOUtSideGC = 0; 
    ICPrepare(NAg_T2,A,AgICprop,ICMemExpMultiplier*ICTypeNums,NBCR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Amongst the cells that were left out and did not enter the GC we
    %%%% select the memory cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     res = SplitMemoryNaive(OutSideCellsArrSp,AgDepClones,Nsp);
    res = SplitMemoryNaive2(OutSideCellsArrSp,AgDepClones,Nsp);
    CellsArrMemorySp = res{1};
    NspArrMemory = res{2};
    OutSideCellsArrSp = res{3};
    
            % Decay of Ag independent memory B cells
    CellArrMemAgInDep = res{4}; % Ag independent memory B cell will decay
    NspArrMemAgInDep = res{5};
    displayCellNum(CellArrMemAgInDep,'Number of Ag independent memory B cells FC');
    res = ExpansionDecay(CellArrMemAgInDep,NspArrMemAgInDep,cellNumDynGrowth,C,dt,TExpansionDecayAgID,Nsp);
    displayCellNum(res{1},'Number of Ag independent memory B cells FC after decay');
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,res{1},Nsp);    
    
    displayCellNum(CellsArrMemorySp,'Number of Ag specific memory B cells FC');
    % Initial memory B cell proliferation
    Res = GrownMemBcells(dt,TMemoryGrowth,NspArrMemory,C,CellsArrSpMemoryExit,CellsArrMemorySp,Nsp,ExitProbOUtSideGC,cellNumDynComp,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt,AgThExp);
    CellsArrMemorySpProf=Res{3};
    displayCellNum(CellsArrMemorySpProf,'Number of Ag specific memory B cells FC following proliferation');
    NspArrMemoryProf = Res{12};
    
    % Merge the B cells created in the GC with those outside
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSpExit,Nsp); 
    displayCellNum(OutSideCellsArrSp,'FC: Naive and mem created in the GC. Before expansion');

    % Merge the proliferated B cell population with the population outside
%     OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    OutSideCellsArrSpAllCells = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProf,Nsp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    %%%%%%%%
    
    NN = 0;
    if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
        for n=1:Nsp
            if(~isempty(OutSideCellsArrSpAllCells{n}))
                idx = find([OutSideCellsArrSpAllCells{n}.Mature]);
                if(length(idx))
                    NN = NN + length(idx);
                    EAgAbExitPop{4,n} = reshape(GetAffField(OutSideCellsArrSpAllCells{n}(idx),'EAbAg'),NAnt+1,length(OutSideCellsArrSpAllCells{n}(idx)))';
                    AgCapExitPop{4,n} = [OutSideCellsArrSpAllCells{n}(idx).AgCap];
%                     MutDist{4,n} = [OutSideCellsArrSpAllCells{n}.ParentNum];
                    BCRsCells = [OutSideCellsArrSpAllCells{n}.BCR_s];
                    MutDist{4,n} = [BCRsCells.ParentNum];
                end
            end
        end
    end
    
    displayCellNum(OutSideCellsArrSpAllCells,'Total Flu challenge. After expansion');
    displayCellNum(CellsArrSpExit,'Exit Flu challenge')
    disp([num2str(NN),' Mature B cells  Flu challenge']);
    displayCellNum(CellsArrSp,'Comp Flu challenge')

    %%%%% The expanded populaiton decayed %%%%%%%   
%     TExpansionDecay = TMemoryGrowth;
    res = ExpansionDecay(CellsArrMemorySpProf,NspArrMemoryProf,cellNumDynGrowth,C,dt,TExpansionDecay,Nsp);
%     res = RemoveExpandedMem(CellsArrMemorySpProf,Nsp)
    
    CellsArrMemorySpProfLeft = res{1};
    displayCellNum(CellsArrMemorySpProfLeft,'Expansion decay memory B cells left FC');
    
    OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrMemorySpProfLeft,Nsp);
    displayCellNum(OutSideCellsArrSp,'Total FC. After expansion decay');
    
    
    end

    %%%%%%%%%%%%%%
    
    InPopNotSeed1 = [];
    InExtNotSeed1 = [];
    InPopNotSeed2 = [];
    InExtNotSeed2 = [];
    InPopNotSeed3 = [];
    InPopNotSeed4 = [];
    C = 1;
    Pop = find(cellNumDynComp(C,:));
    SeedArr1 = SeedingCells{1}(:,1)';
        
    InPopNotSeed1 = [InPopNotSeed1 ; length(setdiff(Pop,SeedArr1))];

   
    C = T1/dt +2;
    
    Pop = find(cellNumDynComp(C,:));
    SeedArr2 = SeedingCells{2}(:,1)';
        
    InPopNotSeed2 = [InPopNotSeed2 ; length(setdiff(Pop,SeedArr2))];


    C = T1/dt+T2/dt+3;
    Pop = find(cellNumDynComp(C,:));
    SeedArr3 = SeedingCells{3}(:,1)';
    
    InPopNotSeed3 = [InPopNotSeed3 ; length(setdiff(Pop,SeedArr3))];
    
    C = T1/dt+T2/dt+T3/dt+4;
    Pop = find(cellNumDynComp(C,:));
    if(length(SeedingCells{4})==0)
        SeedArr4 = [];
    else
        SeedArr4 = SeedingCells{4}(:,1)';
    end
    InPopNotSeed4 = [InPopNotSeed4 ; length(setdiff(Pop,SeedArr4))];
    
    find(InPopNotSeed1)
    find(InPopNotSeed2)
    find(InPopNotSeed3)
    find(InPopNotSeed4)
    
    if( find(InPopNotSeed1) | find(InPopNotSeed2) | find(InPopNotSeed3) |find(InPopNotSeed4))
        InPopNotSeed1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% End of the simulation %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    lambdaSim{k} = lambdaPop;
    AgCapPopSim{k} = AgCapPop;
    AgCapMeanSim{k} = AgCapMean;
    
%     AgCapMeanEpitopesSim{k} = AgCapMeanEpitopes;
%     AgCapEpitopeCell{k,3} = AgCapEpitopeT;
%     AgCapSTDSim{k} = cellfun(@single,AgCapSTD,'un',0);
    AgCapSTDSim{k} = AgCapSTD;
    BirthNumSim{k} = BirthNum;
    MutDistSim{k} = MutDist;
    ICArrtSim{k} = ICArrt;
    EAgAbSim{k} = EAgAbPop;
    EAgAbExitSim{k} = EAgAbExitPop;
    AgCapExitSim{k} = AgCapExitPop;
    %     qSim = {};
    q1Sim{k} = q1Pop;
    q2Sim{k} = q2Pop;
%     PopulationTrack{k,2}=single(cellNumDyn);
    PopulationTrack{k,2}=single(cellNumDynGrowth);
    PopulationTrack{k,2}=single(cellNumDynComp);
    
    SeedingCellsSim{k} = SeedingCells; 
    
    OutEngsCellSim{k} = OutEngsRun;

%     if(exitFlag==0)
%         if(length(TArr{k})==0)
%             TArr{k} = T0+t;
%         else
%             Atemp = TArr{k};
%             Atemp = [Atemp , T0+t];
%             TArr{k} = Atemp;
%         end
%     end
    
    for n=1:Nsp
        cellnumSel(k,n) = length(CellsArrSp{n});
    end
    if(FlagDeathStop && exitFlag)
        for n=1:Nsp
            cellnumSel(k,n) = length(CellsArrSp{n});
        end
        continue;
    end
    
    [IC1tArrAvg IC2tArrAvg] = calc_IC_avg(PopulationTrack,ICArrtSim,k);
    %save([flag_str,'_runs',num2str(runNum),'_mem',num2str(memnumber),'_CorrVV',num2str(Corr12),'mu=',num2str(FitAdv),'w=',num2str(2*Multiplier),'seed=',num2str(Agcutoff),'C=000.mat'],'lambdaSim','AgCapMeanSim','AgCapSTDSim','IC1tArrAvg','IC2tArrAvg','EAgAbSim','EAgAbExitSim','AgCapExitSim',...
     %   'q1Sim','q2Sim','AgCapPopSim','PopulationTrack','BCR_s','EAbAg','EAgMem','wcv','ExitProb','ICTypeNums','AgICprop','runNum',...
      %  'Nsp','T0','T1','T2','T3','Nmax','k','BirthNumSim','MutDistSim','SeedingCellsSim','OutEngsCellSim',...
       % 'BindFlag','paramFlag','NAg','ClusterSize','A','CT','id','OnTargetRCF','OffTargetRCF','FitAdv','bOff','bOn','COff','COn','-v7.3');
    %save(['eng',flag_str,'_runs',num2str(runNum),'_mem',num2str(memnumber),'_CorrVV',num2str(Corr12),'mu=',num2str(FitAdv),'w=',num2str(2*Multiplier),'seed=',num2str(Agcutoff),'C=000.mat'],'OutEngsCellSim');

  save([flag_str,'_runs',num2str(runNum),'_mem',num2str(memnumber),'_CorrVV',num2str(Corr12),'mu=',num2str(FitAdv),'w=',num2str(2*Multiplier),'seed=',num2str(Agcutoff),'mut=0.0.mat'],'lambdaSim','AgCapMeanSim','AgCapSTDSim','IC1tArrAvg','IC2tArrAvg','EAgAbSim','EAgAbExitSim','AgCapExitSim',...
        'q1Sim','q2Sim','AgCapPopSim','PopulationTrack','BCR_s','EAbAg','EAgMem','wcv','ExitProb','ICTypeNums','AgICprop','runNum',...
        'Nsp','T0','T1','T2','T3','Nmax','k','BirthNumSim','MutDistSim','SeedingCellsSim','OutEngsCellSim',...
        'BindFlag','paramFlag','NAg','ClusterSize','A','CT','id','OnTargetRCF','OffTargetRCF','FitAdv','bOff','bOn','COff','COn','-v7.3');
    save(['eng',flag_str,'_runs',num2str(runNum),'_mem',num2str(memnumber),'_CorrVV',num2str(Corr12),'mu=',num2str(FitAdv),'w=',num2str(2*Multiplier),'seed=',num2str(Agcutoff),'mut=0.0.mat'],'OutEngsCellSim');

end

clear aadir
clear RunsDataArr
clear countt
