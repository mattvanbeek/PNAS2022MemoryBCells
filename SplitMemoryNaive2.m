function res = SplitMemoryNaive2(OutSideCellsArrSp,AgDepClones,Nsp)

CellsArrSpNotSelected =cell(1,Nsp);

%%%%%%
%%%%%%
%   OutSideCellsArrSp = MergeCells(OutSideCellsArrSp,CellsArrSpExit,Nsp);
%     NN = 0;
%     if((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
%         for n=1:Nsp
%             if(~isempty(OutSideCellsArrSp{n}))
%                 idx = find([OutSideCellsArrSp{n}.Mature]);
%                 if(length(idx))
%                     NN = NN + length(idx);
%                     EAgAbExitPop{2,n} = reshape(GetAffField(OutSideCellsArrSp{n}(idx),'EAbAg'),3,length(OutSideCellsArrSp{n}(idx)))';
%                     AgCapExitPop{2,n} = [OutSideCellsArrSp{n}(idx).AgCap];
%                 end
%             end
%         end
%     end
%     
%%%%%%
%%%%%%
%%%%%%
%%%%%%

ClonesSizes = cellfun('length',OutSideCellsArrSp);

% TotCellNum = sum(ClonesSizes);

% ProbSeed = min(Nsp/TotCellNum,1);

NspArr = [];
% CloneNum = length(OutSideCellsArrSp);

NspArrMemInAg = [];

CellsArrSp =cell(1,Nsp);
MemoryCellsArrSp = cell(1,Nsp);

CellArrMemAgInDep = cell(1,Nsp);

% for n=1:Nsp

% for n=1:CloneNum
for n=1:Nsp
    if(~isempty(OutSideCellsArrSp{n}))

        if( ismember(n,AgDepClones))
            idxMature = find([OutSideCellsArrSp{n}.Mature]);
            idxNaive = find([OutSideCellsArrSp{n}.Mature]==0);
            
            MemoryClones = OutSideCellsArrSp{n}(idxMature);
            NaiveClones = OutSideCellsArrSp{n}(idxNaive);
            
            CellsArrSp{n} = MemoryClones;
            CellsArrSpNotSelected{n} = NaiveClones;
            if(~isempty(MemoryClones))
                NspArr = [NspArr , n];
            end
        else
            
            idxMature = find([OutSideCellsArrSp{n}.Mature]);
            idxNaive = find([OutSideCellsArrSp{n}.Mature]==0);
            MemoryClonesAgInDep = OutSideCellsArrSp{n}(idxMature); % Choose the population of Ag independent memory B cells that will decay
            CellArrMemAgInDep{n} = MemoryClonesAgInDep;
            NaiveClones = OutSideCellsArrSp{n}(idxNaive);
%             CellsArrSpNotSelected{n} = OutSideCellsArrSp{n};
            CellsArrSpNotSelected{n} = NaiveClones;
            
            if(~isempty(MemoryClonesAgInDep))
                NspArrMemInAg = [NspArrMemInAg , n];
            end
            
        end
    end
end

NspArr = unique(NspArr);
NspArr = sort(NspArr);

res = {};
res{1} = CellsArrSp;
res{2} = NspArr;
res{3} = CellsArrSpNotSelected;
res{4} = CellArrMemAgInDep;
res{5} = NspArrMemInAg;

end