function res = SplitMemoryWeak(OutSideCellsArrSp,AgDepClones,Nsp,Agcutoff)

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
            idxViable=[];
            idxNonViable=[];
            for ii=1:length(OutSideCellsArrSp{n})
          
                    if OutSideCellsArrSp{n}(ii).AgCap >=Agcutoff
                    idxViable=[idxViable;ii];
                    else
                    idxNonViable=[idxNonViable;ii];
                end
            end
%             idxMature = find([OutSideCellsArrSp{n}.Mature]);
%             idxNaive = find([OutSideCellsArrSp{n}.Mature]==0);
            
            ViableClones = OutSideCellsArrSp{n}(idxViable);
            NonViableClones = OutSideCellsArrSp{n}(idxNonViable);
            
            CellsArrSp{n} = ViableClones;
            CellsArrSpNotSelected{n} = NonViableClones;
            if(~isempty(ViableClones))
                NspArr = [NspArr , n];
            end
 
            
        end
    end


NspArr = unique(NspArr);
NspArr = sort(NspArr);

res = {};
res{1} = CellsArrSp;
res{2} = NspArr;
res{3} = CellsArrSpNotSelected;


end