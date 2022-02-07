function OutEngs =  CalcEngMatureNaive(CellsArrSp)
AllClones = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,101,102,103,104,105,106,107,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,133,136,137,138,139,141,142,143,144,145,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,190,195,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228];

OnTargetClonesNPSS = [11 13 14 152 167 169 170 171 172 173]; %NPSS
AllClonesNPSS = [1 3 4 5 6 7 8 9 10 11 12 13 14 91 139 141 144 145 147 148 149 150 151 152 153 155 156 157 158 159 160 161 162 163,...
    164 165 166 167 168 169 170 171 172 173 174 175 190 195 209 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228];

OffTargetClonesNPSS = setdiff(AllClonesNPSS,OnTargetClonesNPSS);

OnTargetClonesNPHA = [11 13 14 128 129 152 153 167 169 170 171 172 173]; % NPHA

OffTargetClonesNPHA = setdiff(AllClones,AllClonesNPSS);


OutEngs = {};

Cell2Check = CellsArrSp;
ClonesSizes = cellfun('length',Cell2Check);
idxPrime = find(ClonesSizes(1,:));
EngsMemory = [];
EngsNaive = [];
Total = 0;

Origin=[];
Seednum=[];
AllEnergy=[];
for ll=1:length(idxPrime)
    CloneCells = Cell2Check{idxPrime(ll)};
        idxNaive = find([CloneCells.Mature]==0);
        idxMature = find([CloneCells.Mature]==1);
        Total = Total + length(idxNaive) + length(idxMature);
            for kk=1:length(CloneCells);
            Origin=[Origin; CloneCells(kk).BCR_s.Origin];
            Seednum=[Seednum;CloneCells(kk).BCR_s.Seednum];
            AllEnergy=[AllEnergy;CloneCells(kk).BCR_s.aff.EAbAg(2:3)];
            end
        
        if(length(idxNaive))
            for kk=1:length(idxNaive)
            EngsNaive = [EngsNaive ; CloneCells(idxNaive(kk)).BCR_s.aff.EAbAg(2), CloneCells(idxNaive(kk)).BCR_s.aff.EAbAg(3)];
            end
        end 
        
        if(length(idxMature))
            for kk=1:length(idxMature)
            EngsMemory = [EngsMemory ; CloneCells(idxMature(kk)).BCR_s.aff.EAbAg(2), CloneCells(idxMature(kk)).BCR_s.aff.EAbAg(3)];
            end
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

EOn = [];
for pp = 1:length(OnTargetClonesNPSS)
    n=OnTargetClonesNPSS(pp);
    HasCells = length([Cell2Check{n}]);
    if(length([Cell2Check{n}]))
        idx = [1:length([Cell2Check{n}])];
        Eclone = reshape(GetAffField(Cell2Check{n}(idx),'EAbAg'),4,length(Cell2Check{n}(idx)))';
        EOn = [EOn ; Eclone(:,2), Eclone(:,3),  Eclone(:,4)];
    end
end

EOff = [];
for pp = 1:length(OffTargetClonesNPSS)
    n=OffTargetClonesNPSS(pp);
    HasCells = length([Cell2Check{n}]);
    if(length([Cell2Check{n}]))
        idx = [1:length([Cell2Check{n}])];
        Eclone = reshape(GetAffField(Cell2Check{n}(idx),'EAbAg'),4,length(Cell2Check{n}(idx)))';
        EOff = [EOff ; Eclone(:,2), Eclone(:,3),Eclone(:,4)];
    end
end

E.on = EOn;
E.off = EOff;

%%%%%%%%%%%%%%
%%%%%%%%%%%%%%


AgCapOn = [];
for pp = 1:length(OnTargetClonesNPSS)
    n=OnTargetClonesNPSS(pp);
    for kk=1:length(Cell2Check{n})
        AgCapOn = [AgCapOn ; Cell2Check{n}(kk).AgCap];
    end
end

AgCapOff = [];
for pp = 1:length(OffTargetClonesNPSS)
    n=OffTargetClonesNPSS(pp);
    for kk=1:length(Cell2Check{n})
        AgCapOff = [AgCapOff ; Cell2Check{n}(kk).AgCap];
    end
end

AgCap.on = AgCapOn;
AgCap.off = AgCapOff;

%%%%%%%%%%%
%%%%%%%%%%%

OutEngs{1} = EngsNaive;
OutEngs{2} = EngsMemory;
OutEngs{3} = E;
OutEngs{4} = AgCap;
OutEngs{5} = [Seednum ,AllEnergy, Origin];
end

