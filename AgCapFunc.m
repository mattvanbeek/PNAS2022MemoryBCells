function [AgCapMeanSimAvg AgCapSTDSimAvg] = AgCapFunc(RunsData)
N = RunsData.runNum;

AgCapMeanSimAvg = zeros(161,1);
AgCapSTDSimAvg = zeros(161,1);

% AgCapMeanSimAvg = RunsData.AgCapMeanSim{1};
% AgCapSTDSimAvg = RunsData.AgCapSTDSim{1};
Count=0;
for i=1:N
    if(length(RunsData.AgCapMeanSim)>=i)
        if( length(RunsData.AgCapMeanSim{i})==161)
            AgCapMeanSimAvg = AgCapMeanSimAvg + RunsData.AgCapMeanSim{i};
            AgCapSTDSimAvg = AgCapSTDSimAvg + RunsData.AgCapSTDSim{i};
            Count = Count+1;
        end
    else
        i
    end
end
AgCapMeanSimAvg = AgCapMeanSimAvg/Count;
AgCapSTDSimAvg = AgCapSTDSimAvg/Count;
end