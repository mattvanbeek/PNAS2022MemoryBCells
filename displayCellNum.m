function displayCellNum(CellArr,str)
Nsp = length(CellArr);
Ntot = 0;
for n=1:Nsp
    Ntot = Ntot + length(CellArr{n});
end
str2print = [str,' ',num2str(Ntot)];
disp(str2print)
end