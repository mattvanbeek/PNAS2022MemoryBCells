function E = GetCrossCell(EAgAbPopt)
E = [];
    for i=1:length(EAgAbPopt)
        E = [E ; EAgAbPopt{i}];
    end
end