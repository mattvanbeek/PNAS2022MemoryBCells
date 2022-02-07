function SeedingCellsArr = SeedingCellsCheck(CellsArrSp,NspArr)
SeedingCellsArr = zeros(length(NspArr),3);
count = 1;
for i=1:length(NspArr)
    n = NspArr(i);
% for n=NspArr
    SeedingCellsArr(count,1) = n;
    SeedingCellsArr(count,2) = length(find([CellsArrSp{n}.Mature]==0));
    SeedingCellsArr(count,3) = length(find([CellsArrSp{n}.Mature]==1));
    count = count + 1;
end
end