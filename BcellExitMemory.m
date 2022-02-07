% function [CellsArr CellsArrExit BirthNum] = BcellExit(CellsArr,lambdaMeanAll,MuFlag,Ntot,qProlCells,ExitProb)
function [CellsArr , CellsArrExit] = BcellExitMemory(CellsArr,ExitProb,MeanAgCapPopulation)
% Exit does not depend on the amount of captured Ag

CellsArrExit = [];
ExitIdx = [];

for i=1:length(CellsArr)
    % Exit Probability does not depend on amount of captured Ag

    if (rand<ExitProb)
        ExitIdx = [ExitIdx ; i];

    end

    
    % Exit Probability depends on amount of captured Ag
%     if(CellsArr(i).AgCap >= MeanAgCapPopulation)
%         if (rand<ExitProb)
%             ExitIdx = [ExitIdx ; i];
%         end
%     else
%         if (rand<ExitProb/2)
%             ExitIdx = [ExitIdx ; i];
%         end
%     end
    
end

idx = [1:1:length(CellsArr)]';

CellsArrExit = CellsArr(ExitIdx);
idxSur = setdiff(idx,ExitIdx);

CellsArr = CellsArr(idxSur);


end