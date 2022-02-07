function ret = GetAffField(BcellArr,f)
a = [BcellArr.BCR_s];
b = [a.aff];
ret = [b.(f)];
end