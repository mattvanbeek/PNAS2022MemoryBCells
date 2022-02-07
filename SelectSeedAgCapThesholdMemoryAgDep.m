% function res = SelectSeedAgCapThMemory(CellsArrSpExit,EAgAbExitPop,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,MedianAgCapThreshold)
function res = SelectSeedAgCapThesholdMemoryAgDep(CellsArrSpExit,Nsp,paramFlag,EnterProbBelowMed,EnterProbAboveMed,AgThresholdGCSeed,AgDepClones)

EnergySeedClones = 0;
% if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
%     for n=1:Nsp
%         if(~isempty(CellsArrSpExit{n}))
%             EAgAbExitPop{1,n} = reshape(GetAffField(CellsArrSpExit{n},'EAbAg'),3,length(CellsArrSpExit{n}))';
%             AgCapExitPop{1,n} = [CellsArrSpExit{n}.AgCap];
%         end
%     end
% end


%%%%%%%%%


ExitCellEnergy = [];
for n=1:length(CellsArrSpExit)
    n

    EnergypExitClone = CellsArrSpExit{n};
    if(~isempty(EnergypExitClone))
        ExitCellEnergy = [ExitCellEnergy ; reshape(GetAffField(CellsArrSpExit{n},'EAbAg'),4,length(CellsArrSpExit{n}))'];
        
    end
end
ExitCellEnergy = ExitCellEnergy(:,2);

%%%%%%%%%

ExitCellAgCap = [];
for n=1:length(CellsArrSpExit)
    AgCapExitClone = CellsArrSpExit{n};
    if(~isempty(AgCapExitClone))
        ExitCellAgCap = [ExitCellAgCap  [AgCapExitClone.AgCap]];
    end
end

%%%%%%%%%%%%%%

MedianEng = median(ExitCellEnergy);
MedianAgCap = median(ExitCellAgCap);

%%%%%%%%%%%%%%


%%%%%%%%%%%%%%

NspArr = [];
CellsArrSp =cell(1,Nsp);
OutSideCellsArrSp =cell(1,Nsp);
for n=1:Nsp
    if(~isempty(CellsArrSpExit{n}))
        InputClone = [];LeftClones = [];
        for j=1:length(CellsArrSpExit{n})
            
            if(EnergySeedClones)
                % Select based on the energy
                EAbAgCell = CellsArrSpExit{n}(j).BCR_s.aff.EAbAg(2);
                FlagTake = 0;
                if(EAbAgCell<=MedianEng)
                    if(rand<EnterProbBelowMed)
                        FlagTake = 1;
                    end
                else
                    if(rand<EnterProbAboveMed)
                        FlagTake = 1;
                    end
                end
            else
                % Select based on the amount of captured Ag
                AgCapCell = CellsArrSpExit{n}(j).AgCap;
                FlagTake = 0;
                if(AgCapCell<AgThresholdGCSeed)
                    if(rand<EnterProbBelowMed)
                        FlagTake = 1;
                    end
                else
                    if(rand<EnterProbAboveMed)
                        FlagTake = 1;
                    end
                end
                
            end
            if( ismember(n,AgDepClones))
                if(FlagTake)
                    InputClone = [InputClone ; CellsArrSpExit{n}(j)];
                else
                    LeftClones = [LeftClones ; CellsArrSpExit{n}(j)];
                end
            else
                LeftClones = [LeftClones ; CellsArrSpExit{n}(j)];
            end
        end
        CellsArrSp{n} = InputClone;
        OutSideCellsArrSp{n} = LeftClones;
        if(~isempty(InputClone))
            NspArr = [NspArr , n];
        end
    end
end

NspArr = unique(NspArr);
NspArr = sort(NspArr);

res = {};
res{1} = CellsArrSp;
res{2} = NspArr;
res{3} = OutSideCellsArrSp; % These are the cells that were left out and did not enter the GC reaction

end