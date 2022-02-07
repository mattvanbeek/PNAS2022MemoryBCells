function Res = GrownMemBcells(dt,T,NspArr,C,CellsArrSpExit,CellsArrSp,Nsp,ExitProb,cellNumDyn,AgCapMean,AgCapSTD,lambdaPop,AgCapPop,ICArrt,AgThExp)
global ICArr;
MuFlag=2;

% AgThExp = 100;

if(length(NspArr)==0)
    
    Res = {};
    Res{1}=C;
    Res{2}=cell(1,Nsp);
    Res{3}=CellsArrSp;
    Res{4}=Nsp;
    Res{5}=ExitProb;
    Res{6}=cellNumDyn;
    Res{7}=[];
    Res{8}=[];
    Res{9}=cell(1,Nsp);
    Res{10}=cell(1,Nsp);
    Res{11}={};
    Res{12}=NspArr;
    return;
end


for t=0:dt:T
        
        Ntot = 0;
        
        % Accumulate Exit B cells
%         CellsArrSpExit = AccumulateExitBCs(CellsArrSpExit,CellsArrSp,Nsp,ExitProb,0);
        CellsArrSpExit = AccumulateExitBCsMark(CellsArrSpExit,CellsArrSp,Nsp,ExitProb,0);
        for n=1:Nsp
            Ntot = Ntot + length(CellsArrSp{n});
        end
        
        idxSps = randperm(length(NspArr)); % randomize the order of cells extracting Ags
        for n=NspArr(idxSps)
            if(length(CellsArrSp{n})==0)
                NspArr = NspArr(find(NspArr ~= n));
                %                 exitFlag = 1;
            end
        end
        
        [lambdaMeanAll AgCapMeanAll CTMeanAll] = LambdaMeanFuncAgCapHABoost2(CellsArrSp,NspArr,MuFlag,AgThExp);
        
        exitFlag = 0;
        qProlCell_n = [];
        idxSps = randperm(length(NspArr)); % randomize the order of cells extracting Ags
        for n=NspArr(idxSps)
            
            [CellsArrSp{n} qProlCell BN MeanCapAg] = BD_step(CellsArrSp{n},lambdaMeanAll,AgCapMeanAll,CTMeanAll,MuFlag,Ntot);
%             AgCapEpitopeT(C,n) = MeanCapAg;
            BirthNum{C,n} = single(BN);
            MutDist{C,n} = single([CellsArrSp{n}.ParentNum]);
            if(length(CellsArrSp{n})==0)
                NspArr = NspArr(find(NspArr ~= n));
                exitFlag = 1;
            end
            qProlCell_n = [qProlCell_n ; qProlCell];
        end
        
        
        % q(1) : the kon of the conserved site ; q(2) : the kon to the variable site of Ag1 ; q(3) : the kon to the variable site of Ag2
%         qProlCells_t{Countq} = qProlCell_n;  % qa0 = obj.aff.wcv(1)*obj.q(1)+obj.aff.wcv(2)*obj.q(2); qb0 = obj.aff.wcv(1)*obj.q(1)+obj.aff.wcv(2)*obj.q(3);
%         Countq = Countq+1;
        %         qProlCells_t{C} = qProlCells;
        %%%%%
        
        
        % paramFlag=1 Mutate EAbAg
        % paramFlag=2 Mutate q1
        % paramFlag=3 Mutate q2
        
        AgCapMean(C,1) = AgCapMeanAll;
        AgCapSTD(C,1) = std(AgCapDist(CellsArrSp,NspArr));
        
        for n=NspArr
            lambdaPop{C,n} = single([CellsArrSp{n}.lambda]);
            AgCapPop{C,n} = single([CellsArrSp{n}.AgCap]);
%             
%             if(paramFlag==1)
%                 EAgAbPop{C,n} = reshape(GetAffField(CellsArrSp{n},'EAbAg'),3,length(CellsArrSp{n}))';
%             elseif(paramFlag==2)
%                 q1Pop{C,n} = reshape(GetAffField(CellsArrSp{n},'q1'),3,length(CellsArrSp{n}))';
%             elseif(paramFlag==3)
%                 q2Pop{C,n} = GetAffField(CellsArrSp{n},'q2');
%                 Agcap_n = [CellsArrSp{n}.AgCap];
%             elseif(paramFlag==4)
%                 EAgAbPop{C,n} = GetAffField(CellsArrSp{n},'EAbAg');
%                 q1Pop{C,n} = GetAffField(CellsArrSp{n},'q1');
%             elseif((paramFlag==5) || (paramFlag==7) || (paramFlag==8))
%                 EAgAbPop{C,n} = reshape(GetAffField(CellsArrSp{n},'EAbAg'),3,length(CellsArrSp{n}))';
%             elseif(paramFlag==6)
%                 q1Pop{C,n} = reshape(GetAffField(CellsArrSp{n},'q1'),3,length(CellsArrSp{n}))';
%             end
        end
        
        for n=NspArr
            cellNumDyn(C,n) = length(CellsArrSp{n});
        end
        
        ICArrt{2,C} = ICArr; %??????
        
        %%%%%
        C = C+1;
        
        if(exitFlag)
%             if(length(TArr{k})==0)
%                 TArr{k} = T0+t;
%                 if(FlagDeathStop)
%                     break;
%                 end
%             else
%                 Atemp = TArr{k};
%                 Atemp = [Atemp , T0+t];
%                 TArr{k} = Atemp;
%             end
            if(length(NspArr)==0)
                break;
            end
        end
end


Res = {};
Res{1}=C;
Res{2}=CellsArrSpExit;
Res{3}=CellsArrSp;
Res{4}=Nsp;
Res{5}=ExitProb;
Res{6}=cellNumDyn;
Res{7}=AgCapMean;
Res{8}=AgCapSTD;
Res{9}=lambdaPop;
Res{10}=AgCapPop;
Res{11}=ICArrt;
Res{12}=NspArr;

end