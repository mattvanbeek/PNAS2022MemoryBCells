function [lambdaMeanAll AgCapMeanAll CTMeanAll] = LambdaMeanFuncAgCapHABoost2(CellsArrSp,NspArr,MuFlag,AgThExp)

% lambda_all =[];
AgCapCloneAll =[];
for n=NspArr
    AgCapClone = [CellsArrSp{n}.AgCap];
%     idx = ones(1,length(AgCapClone));
%     lambdanEff = idxMul.*lambdan;
%     lambda_all = [ lambda_all , lambdanEff];
    AgCapCloneAll = [AgCapCloneAll , AgCapClone];
end


if(MuFlag==1)
idxlargezero = find(AgCapCloneAll); %in the average we take only B cells the captured at least one Ag molecule
% lambdaMeanAll = mean(lambda_all);
if(length(idxlargezero))
%     lambdaMeanAll = mean(lambda_all(idxlargezero));
    AgCapMeanAll = mean(AgCapCloneAll(idxlargezero));
else
    AgCapMeanAll=0;
end
elseif(MuFlag==2)
    idxlargezero = find(AgCapCloneAll>=AgThExp); %in the average we take only B cells the captured at least one Ag molecule
% lambdaMeanAll = mean(lambda_all);
if(length(idxlargezero))
%     lambdaMeanAll = mean(lambda_all(idxlargezero));
    AgCapMeanAll = mean(AgCapCloneAll(idxlargezero));
else
    AgCapMeanAll=0;
end
end

lambda = [];
CTCells = [];
for n=NspArr
    AgCapClone = [CellsArrSp{n}.AgCap];
    bArr = [];
    CTArr = [];
    
    for i=1:length(CellsArrSp{n})
        bArr(i) = CellsArrSp{n}(i).BCR_s.b;
        CTArr(i) = CellsArrSp{n}(i).BCR_s.aff.CT;
    end
%     lambdaBaseArr = [CellsArrSp{n}.lambdaBase];
    CTCells = [CTCells, CTArr];
%     lambda = [lambda , bArr.*lambdaBaseArr.*(CTArr + AgCapClone)./(CTArr + AgCapMeanAll)];
%     [CellsArrSp{4}.lambdaBase].*
%     obj.lambdaBase*( (CT + obj.AgCap)/(CT + lambdaMeanAll));
end
CTMeanAll = mean(CTCells);
for n=NspArr
    AgCapClone = [CellsArrSp{n}.AgCap];
    bArr = [];
    CTArr = [];
    
    for i=1:length(CellsArrSp{n})
        bArr(i) = CellsArrSp{n}(i).BCR_s.b;
        CTArr(i) = CellsArrSp{n}(i).BCR_s.aff.CT;
    end
    lambdaBaseArr = [CellsArrSp{n}.lambdaBase];
%     CTCells = [CTCells, CTArr];
    lambda = [lambda , bArr.*lambdaBaseArr.*(CTArr + AgCapClone)./(CTMeanAll + AgCapMeanAll)];
%     [CellsArrSp{4}.lambdaBase].*
%     obj.lambdaBase*( (CT + obj.AgCap)/(CT + lambdaMeanAll));
end


lambdaMeanAll = mean(lambda);


end