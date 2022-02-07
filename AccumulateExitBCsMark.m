function CellsArrSpExit = AccumulateExitBCsMark(CellsArrSpExit,CellsArrSp,Nsp,ExitProb,lambdaMeanAll)

for n=1:Nsp
    [CellsArrSp{n} , CellsArrExit] = BcellExitMemory(CellsArrSp{n},ExitProb,lambdaMeanAll);
    for j=1:(length(CellsArrExit))
        CellsArrExit(j).Mature = 1;
    end
    if(length(CellsArrSpExit) <n)
        CellsArrSpExit{n} = CellsArrExit;
    else
        temp = CellsArrSpExit{n};
        temp = [temp ; CellsArrExit];
        CellsArrSpExit{n} = temp;
    end
end

end