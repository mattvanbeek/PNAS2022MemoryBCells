function res = GrowthPhase(CellsArrSp,NspArr,cellNumDyn,C,dt,T0,Nsp)
global ICArr;
% global ICArr;
% global ICArr;

if(length(NspArr)==0)
    exitFlag = 1;
    res={};
    res{1} = CellsArrSp;
    res{2} = NspArr;
    
%     cellNumDyn = [];
%     for i=1:length(NspArr)
%         n=NspArr(i);
%         cellNumDyn(C,n) = length(CellsArrSp{n});
%     end
    res{3} = cellNumDyn;
    res{4} = C;
    return;
end




% CT is the same for all
CTMeanAll = 0;

%%% G1 %%%%%%%%

MuFlag = 0;
for t=0:dt:T0
    qProlCells = {};
    t;
    Ntot = 0;
    for n=1:Nsp
        Ntot = Ntot + length(CellsArrSp{n});
    end
    exitFlag = 0;
    
    lambdaMeanAll = LambdaMeanFunc(CellsArrSp,NspArr);
    
    for n=NspArr
%         [CellsArrSp{n} qProlCell MeanCapAg] = BD_step(CellsArrSp{n},lambdaMeanAll,MuFlag,Ntot);
        [CellsArrSp{n} qProlCell MeanCapAg] = BD_step(CellsArrSp{n},lambdaMeanAll,CTMeanAll,0,MuFlag,Ntot);
        
        if(length(CellsArrSp{n})==0)
            NspArr = NspArr(find(NspArr ~= n));
%             if(FlagDeathStop)
%                 exitFlag = 1;
%                 break;
%             else
%                 Atemp = TArr{k};
%                 Atemp = [Atemp , t];
%                 TArr{k} = Atemp;
%             end
            if(length(NspArr)==0)
                exitFlag = 1;
                break;
            end
        end
    end
    
    for n=NspArr
        cellNumDyn(C,n) = length(CellsArrSp{n});
    end
    
    ICArrt{1,C} = ICArr;
    
    C = C + 1;
    
    if(exitFlag)
        break;
    end
end

res={};
res{1} = CellsArrSp;
res{2} = NspArr;
res{3} = cellNumDyn;
res{4} = C;

% PopulationTrack{k,1}=single(cellNumDyn);

% for n=NspArr
%     cellnumGrowth(k,n) = length(CellsArrSp{n});
% end

%%%%%%%%%%%
%% %%%%%%%%%

% if(length(NspArr)==0)
%     continue;
% end

%%%%%%%%%%%
%%%%%%%%%%%
end