function AgDepClones = AgDepClonesList(Ag)
% AllClones = [1:1:228];
% OnTargetClonesNPSS = [11 13 14 152 167 169 170 171 172 173]; %NPSS
% AllClonesNPSS = [1 3 4 5 6 7 8 9 10 11 12 13 14 91 139 141 144 145 147 148 149 150 151 152 153 155 156 157 158 159 160 161 162 163,...
%     164 165 166 167 168 169 170 171 172 173 174 175 190 195 209 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228];
% OffTargetClonesNPSS = setdiff(AllClonesNPSS,OnTargetClonesNPSS);
% OnTargetClonesNPHA = [11 13 14 128 129 152 153 167 169 170 171 172 173]; % NPHA
% OffTargetClonesNPHA = setdiff([1:1:228],AllClonesNPSS);

AllClones = [1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,101,102,103,104,105,106,107,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,133,136,137,138,139,141,142,143,144,145,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,190,195,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228];
OnTargetClonesNPSS = [11 13 14 152 167 169 170 171 172 173]; %NPSS
AllClonesNPSS = [1 3 4 5 6 7 8 9 10 11 12 13 14 91 139 141 144 145 147 148 149 150 151 152 153 155 156 157 158 159 160 161 162 163,...
    164 165 166 167 168 169 170 171 172 173 174 175 190 195 209 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228];
OffTargetClonesNPSS = setdiff(AllClonesNPSS,OnTargetClonesNPSS);
OnTargetClonesNPHA = [11 13 14 128 129 152 153 167 169 170 171 172 173]; % NPHA
OffTargetClonesNPHA = setdiff(AllClones,AllClonesNPSS);


NP1s = '1HA';NP2s = 'NPSS';NP3s = 'NPSS';NP4s = '40VirusHA';
AgDepClones = [];

switch Ag
    case '1HA'
        AgDepClones = AllClones;
    case '40VirusHA'
        AgDepClones = AllClones;
    case 'NPHA'
        AgDepClones = AllClones;
    case 'NPSS'
        AgDepClones = AllClonesNPSS;
end

end