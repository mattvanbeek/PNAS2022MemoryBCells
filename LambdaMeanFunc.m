function lambdaMeanAll = LambdaMeanFunc(CellsArrSp,NspArr)
lambda_all =[];
alpha=2;
for n=NspArr
    idxNL = [CellsArrSp{n}.Flagln];
    lambdan = [CellsArrSp{n}.lambda];
%     lambdan = [CellsArrSp{n}.lambdaPersonal];
    idxMul = ones(1,length(lambdan));
    idxMul(find(idxNL))=alpha;
    lambdanEff = idxMul.*lambdan;
    lambda_all = [ lambda_all , lambdanEff];
% 
%     idxNL = [CellsArrSp{n}.Flagln];
%     lambdan = [CellsArrSp{n}.AgCap];
%     idxMul = ones(1,length(lambdan));
%     idxMul(find(idxNL))=alpha;
%     lambdanEff = idxMul.*lambdan;
%     lambda_all = [ lambda_all , lambdanEff];

end
lambdaMeanAll = mean(lambda_all);
end