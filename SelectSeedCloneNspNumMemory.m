function res = SelectSeedCloneNspNumMemory(CellsArrSpExit,Nsp)
% Nsp is the initial number of clone we now select a subset of size Nsp
% from CellsArrSpExit

% Nsp now is the initial number cell that will seed the GC. We now select a subset of size Nsp
% from CellsArrSpExit

ClonesSizes = cellfun('length',CellsArrSpExit);

TotCellNum = sum(ClonesSizes);

ProbSeed = min(Nsp/TotCellNum,1);

NspArr = [];
CloneNum = length(CellsArrSpExit);

CellsArrSp =cell(1,CloneNum);
MemoryCellsArrSp = cell(1,CloneNum);
% for n=1:Nsp
for n=1:CloneNum
    if(~isempty(CellsArrSpExit{n}))
        InputClone = [];
        LeftClones = [];
        for j=1:length(CellsArrSpExit{n})
            
            FlagTake = 0;
            
            if(rand<ProbSeed)
                FlagTake = 1;
            end
            if(FlagTake)
                InputClone = [InputClone ; CellsArrSpExit{n}(j)];
            else
                LeftClones = [LeftClones ; CellsArrSpExit{n}(j)];
            end
            
        end
        CellsArrSp{n} = InputClone;
        MemoryCellsArrSp{n} = LeftClones;
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
res{3} = MemoryCellsArrSp;
end