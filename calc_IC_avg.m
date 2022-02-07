function [IC1tArrAvg IC2tArrAvg] = calc_IC_avg(PopulationTrack,ICArrtSim,runNum)

IC1tArrp = {};
IC2tArrp = {};
endoftimesArr = [];

for p=1:runNum
    p;
    Pop = PopulationTrack{p,2};
    ICArrt = ICArrtSim{p};
%     endoftimes = size(Pop,1);
    endoftimes = size(ICArrt,2);
    endoftimesArr(p) = endoftimes;
    IC1tArr = [];
    IC2tArr = [];
    for ic_t = 1:endoftimes
%         if(ic_t==26)
%             ic_t
%         end
        ic_t;
        ict1 = [];
        if(~isempty(ICArrt{2,ic_t}))
%             for l=1:size(ICArrt{2,ic_t}{1},2)
%                 l;
%                 for m=1:size(ICArrt{2,ic_t}{1}{l},2)
%                     m;
%                     ict1 = [ict1 ; ICArrt{2,ic_t}{1}{l}{m}];
%                 end
%             end
%             ICt1 = mean(ict1,1);
            
            AA = ICArrt{2,ic_t}{1}';
            out=cell2mat(cellfun(@(x) cell2mat(x),AA,'un',0));
            ICnum = length(AA);
            BCRNum = length(AA{1});
            outreshape = reshape(out',2,ICnum*BCRNum)';
%             ICt1_mean = mean(outreshape);
            ICt1 = mean(outreshape);
%             if(length(find(ICt1_mean-ICt1)))
%                 ICt1;
%             end
                        
%             ict2 = [];
%             for l=1:size(ICArrt{2,ic_t}{2},2)
%                 l;
%                 for m=1:size(ICArrt{2,ic_t}{2}{l},2)
%                     m;
%                     ict2 = [ict2 ; ICArrt{2,ic_t}{2}{l}{m}];
%                 end
%             end
%             ICt2 = mean(ict2,1);
            
            AA = ICArrt{2,ic_t}{2}';
            out=cell2mat(cellfun(@(x) cell2mat(x),AA,'un',0));
            ICnum = length(AA);
            BCRNum = length(AA{1});
            outreshape = reshape(out',2,ICnum*BCRNum)';
%             ICt2_mean = mean(outreshape);
            ICt2 = mean(outreshape);
%             if(length(find(ICt2_mean-ICt2)))
%                 ICt2;
%             end
            
        else
            ICt1=[0 0];
            ICt2=[0 0];
        end
        
        %     ICt1 = ICArrt{2,ic_t}{1}{1}{10};
        %     ICt2 = ICArrt{2,ic_t}{2}{1}{10};
        IC1tArr = [IC1tArr ; ICt1];
        IC2tArr = [IC2tArr ; ICt2];
    end
    IC1tArrp{p} = IC1tArr;
    IC2tArrp{p} = IC2tArr;
end

MaxTime = max(endoftimesArr);
IC1tArrAvg = zeros(MaxTime,2);
IC2tArrAvg = zeros(MaxTime,2);

if( (length(IC1tArrAvg)==0) | (length(IC2tArrAvg)==0))
    IC1tArrAvg=[];
    IC2tArrAvg=[];
    return;
end
IC1tCnt = zeros(MaxTime,1);
IC2tCnt = zeros(MaxTime,1);
for p=1:runNum
    IC1tArrAvg(1:endoftimesArr(p),:) = IC1tArrAvg(1:endoftimesArr(p),:) + IC1tArrp{p};
    IC2tArrAvg(1:endoftimesArr(p),:) = IC2tArrAvg(1:endoftimesArr(p),:) + IC2tArrp{p};
    IC1tCnt(1:endoftimesArr(p)) = IC1tCnt(1:endoftimesArr(p)) + ones(endoftimesArr(p),1);
    IC2tCnt(1:endoftimesArr(p)) = IC2tCnt(1:endoftimesArr(p)) + ones(endoftimesArr(p),1);
end
IC1tArrAvg = IC1tArrAvg./IC1tCnt;
IC2tArrAvg = IC2tArrAvg./IC2tCnt;
end