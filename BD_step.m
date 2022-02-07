function [CellsArr qProlCell BirthNum MeanCapAg] = BD_step(CellsArr,lambdaMeanAll,AgCapMeanAll,CTMeanAll,MuFlag,Ntot,qProlCells)
Deadidx = [];
qProlCell = [];
ThresholdAgDeathCount = 0;
BirthNum = 0;
%         AgThreshold = 10;
AgThreshold = 0;
ZeroAgCaptureCount = 0;
CaptureAttempts = 1;
NumberCells = length(CellsArr);

TotAgCaptured = [];

for i=1:length(CellsArr)
    
%     if(MuFlag == 0)
%         Kid = CellsArr(i).BirthTX();
%     end
%     
%     
    
    if(MuFlag == 0)
        Kid = CellsArr(i).BirthTX();
    else
        if(MuFlag==1 || MuFlag==2)
            CellsArr(i).AgCap = 0;
            for kk=1:CaptureAttempts
%                 i;
%                 CellsArr(i);
%                 CellsArr(i).AffinityClass.CAgCaptureFunc();
%                 CellsArr(i).AgCap;
                CellsArr(i).AgCap = CellsArr(i).AgCap + CellsArr(i).AffinityClass.CAgCaptureFunc();
                TotAgCaptured = [TotAgCaptured ; CellsArr(i).AgCap];
            end
            if(CellsArr(i).AgCap==0)
                ZeroAgCaptureCount = ZeroAgCaptureCount+1;
            end
        end
        if(MuFlag==1)
            Kid = CellsArr(i).BirthNorm(AgCapMeanAll,CTMeanAll);
        elseif(MuFlag==2)
            Kid = CellsArr(i).BirthNormMemoryOutsideGC(AgCapMeanAll,0);
        elseif(MuFlag==3)
            Kid = [];
        end
    end
    
    if(length(Kid))
        qProlCell = [qProlCell ; CellsArr(i).BCR_s.aff.q1 ] ;
        CellsArr(length(CellsArr)+1,1) = Kid;
        BirthNum = BirthNum+1;
    end
    
    if(MuFlag == 0)
        Dflag = CellsArr(i).DeathT();
    elseif(MuFlag == 1)
        Dflag = CellsArr(i).DeathTXn(lambdaMeanAll,Ntot);
        
        if(CellsArr(i).AgCap<AgThreshold)
            Dflag = 1;
            ThresholdAgDeathCount = ThresholdAgDeathCount+1;
        end
    elseif(MuFlag == 2)
        %Dflag = CellsArr(i).DeathT();
%         Dflag = 0;
        Dflag = CellsArr(i).DeathExpansionDecay();
    elseif(MuFlag == 3)
%         Dflag = CellsArr(i).DeathT();
        Dflag = CellsArr(i).DeathExpansionDecay();
    end
    if(Dflag)
        Deadidx = [Deadidx ; i];
    end
end
ZeroAgCaptureCount/NumberCells;

% ThresholdAgDeathCount
if(ThresholdAgDeathCount>1)
    ThresholdAgDeathCount;
end
idx = [1:1:length(CellsArr)]';
idxSur = setdiff(idx,Deadidx);
CellsArr = CellsArr(idxSur);
if(length(TotAgCaptured))
    MeanCapAg = mean(TotAgCaptured);
else
    MeanCapAg=0;
end
% qProlCells{length(qProlCells)+1} = qProlCell;
end
