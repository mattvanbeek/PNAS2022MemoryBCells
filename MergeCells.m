 function CellA = MergeCells(CellA,CellB,Length)
% Merge B into A
% Length is the desired find length of CellA
for n=1:Length
    if(length(CellA) <n)
        CellA{n} = [];
    else
        temp = CellA{n};
        temp = [temp ; CellB{n}];
        CellA{n} = temp;
    end
end

end
