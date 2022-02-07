function ICPrepare(MAg,A,AgICprop,ICTypeNums,NBCR)
% ICArr is a global
global ICArr;
% obj.MAg = MAg;
R = 100;
p = (A/R)^2;

for i=1:2
    Z = sum(AgICprop(i,:));
    MAg_a = round(MAg*AgICprop(i,1)/Z);
    MAg_b = MAg-MAg_a;
    
    IC = {};
    for j=1:ICTypeNums(i)
        for k=1:NBCR
%             NAg_a = binornd(MAg_a,p);
%             NAg_b = binornd(MAg_b,p);
            NAg_a = MAg_a;
            NAg_b = MAg_b;
            AgDens = [NAg_a NAg_b];
            IC{j}{k} = AgDens;
        end
    end
    ICArr{i} = IC;
end
end