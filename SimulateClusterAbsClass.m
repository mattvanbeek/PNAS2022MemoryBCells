classdef SimulateClusterAbsClass < handle
    properties
        NBCR;
        NAg;
        W;
        q;
        tau;
        F;
        EAgMem;
        EAbAg;
        TCstates;
        ACstates;
        StatesNum;
        AgCapPro;
        AgCapProStartingState;
        StartingState;
        StartingStateNames;
        Parr;
        ClusterSize;
    end
    methods
        function obj = SimulateClusterAbsClass(NAg,q,EAgMem,EAbAg,tau,F,ClusterSize)
            if(nargin > 0)
                obj.ClusterSize = ClusterSize;
                obj.NBCR = obj.ClusterSize;
                obj.NAg = NAg;
                
                %                 obj.q=0.2;
                %                 obj.EAbAg = 0.05;
                %                 obj.EAgMem = 1;
                obj.q=q;
                obj.EAbAg = EAbAg;
                              global Virustype
                if Virustype==3
                    obj.EAbAg=[EAbAg(1) EAbAg(2) EAbAg(4) EAbAg(3)];
                    keyboard;
                else
                end
                obj.EAgMem = EAgMem;
                obj.tau=tau;
                obj.F=F;
                
                %                 obj.ConstructWMat();
                %                 obj.CaptureProb();
                %                 obj.StartingStatesFind(obj.TCstates);
                %                 obj.CalcProbDensityNoneq();
                %                 obj.AgCapProStartingState = obj.AgCapPro(obj.StartingState,:);
                
                %                 obj.AbsorbingStatesFind();
                obj.StartingStatesFind(obj.TCstates);
                obj.CalcProbDensityNoneq();
                
            else
            end
        end
        
        
        function StartingState = StartingStatesFind(obj,TCstates)
            NBCR = obj.NBCR;
            
            %             StartingState = [];
            %             StartingStateNames = {};
            %             for i=1:length(TCstates)
            %
            % %                 ClusterState = dec2base(TCstates(i)-1,obj.StatesNum,obj.NBCR);
            % %
            % %                 BoundArms = regexp(ClusterState,'[0-1]');
            % %                 NonBoundAbs = regexp(ClusterState,'[5]');
            %
            %
            %                 ByteCluster = zeros(1,NBCR);
            %                 b=i-1;
            %                 for p=1:NBCR
            %                     ByteCluster(NBCR-p+1) = mod(b,9);
            %                     b = floor(b/9);
            %                 end
            %                 BoundArms = find( (ByteCluster==0) | (ByteCluster==1));
            %                 NonBoundAbs = find(ByteCluster==5);
            %
            %
            %                 if( (length(BoundArms) == 1) & (length(NonBoundAbs) == obj.NBCR-1))
            %                     StartingState = [StartingState , TCstates(i)];
            %                     c = length(StartingStateNames);
            % %                     StartingStateNames{c+1} = ClusterState;
            %                     StartingStateNames{c+1} = ByteCluster;
            %                 end
            %             end
            %             obj.StartingState = StartingState;
            %             obj.StartingStateNames = StartingStateNames;
            
            Allbits = [1:NBCR];
            StartingStateNames1 = {};
            StartingState1 = [];
            ConnectedBits = [0 1];
            count = 1;
            
            for i=1:NBCR
                for j=1:length(ConnectedBits)
                    State = [];
                    State(i) = ConnectedBits(j);
                    idx = setdiff(Allbits,i);
                    State(idx) = 5;
                    StartingStateNames1{count} = State;
                    
                    
                    StateNum = 0;
                    for k=1:NBCR
                        StateNum = StateNum+State(k)*9^(NBCR-k);
                    end
                    StateNum = StateNum+1;
                    StartingState1(count) = StateNum;
                    
                    count = count+1;
                end
            end
            
            obj.StartingState = StartingState1;
            obj.StartingStateNames = StartingStateNames1;
            
        end
        
        function AbsorbingStatesFind(obj)
            
            Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            StatesNum = length(Tstate)+length(Astate);
            obj.StatesNum = StatesNum;
            %             W = sparse(StatesNum^NBCR,StatesNum^NBCR);
            NBCR = obj.NBCR;
            ByteCluster = zeros(1,NBCR);
            for i=0:(StatesNum^NBCR-1)
                %                 ClusterState = dec2base(i,StatesNum,NBCR);
                
                b=i;
                for p=1:NBCR
                    ByteCluster(NBCR-p+1) = mod(b,9);
                    b = floor(b/9);
                end
                
                %                 TransitionIndex = regexp(ClusterState,'[0-4]');
                %                 TransitionNonConIndex = regexp(ClusterState,'[5-7]');
                
                TransitionIndex = find(ByteCluster<=4);
%                 TransitionNonConIndex = find( (ByteCluster>=5) & (ByteCluster<7) )
                TransitionNonConIndex = find( (ByteCluster>=5) & (ByteCluster<=7) );
                
                NumConnected = length(TransitionIndex);
                NumNonConnected = length(TransitionNonConIndex);
                if(NumConnected)
                    TCstates = [TCstates , i+1];
                else
                    ACstates = [ACstates , i+1];
                end
            end
            obj.TCstates = TCstates;
            obj.ACstates = ACstates;
        end
        
        %%%%%
        %%%%%
        
        
        %%%%%
        %%%%%
        
        function CalcProbDensityNoneq(obj);
            
            %             w = obj.w;
            tau = obj.tau;
            n=obj.NAg;
            % We start from 00 and have a chance of getting to 10/01/11
            % before time tau;
            q = obj.q;
            switch n
                case 0
                    obj.Parr(1) = 1;
                    obj.Parr(2) = 0;
                    obj.Parr(3) = 0;
                    obj.Parr(4) = 0;
                otherwise
                    obj.Parr(1) = exp(-2*q*obj.NBCR*obj.NAg*tau);
                    POneArm = 1-obj.Parr(1);
                    p=POneArm/(2*obj.NBCR);
                    for j=1:length(obj.StartingState)
                        obj.Parr(j+1) = p;
                    end
                    %                     obj.Parr(2) = 0.5*(1-exp(-2*q*n*tau));
                    %                     obj.Parr(3) = 0.5*(1-exp(-2*q*n*tau));
                    %                     obj.Parr(4) = 0;
            end
        end
        
        %%%%
        %%%%
        %%%%
        
        function [EndStateName EndStateNum ACstates_AgCap count]  = SimulateClusterAbs(obj,InitialState)
            %             ByteCluster = InitialState;
            %             ss = '0';
            %             for j=1:length(InitialState)
            %                 ss(j) = num2str(InitialState(j));
            %             end
            %             InitialState = ss;
            %             NextStateNum = base2dec(InitialState,obj.StatesNum);
            
            NBCR = obj.NBCR;
            %             NextStateNum = 0;
            %             for i=1:NBCR
            %                 NextStateNum = NextStateNum+InitialState(i)*9^(NBCR-i);
            %             end
            %             NextStateNum = NextStateNum+1;
            
            NextStateName = InitialState;
            count = 0;
            
            while ( isAbsorbing(obj,NextStateName) == 0)
                %             while ( length(find(obj.ACstates == NextStateNum+1)) == 0)
                %                 [NextStateNameArr NextStateNumArr NextStateProbArr] = obj.TransitionSimulate(NextStateName,NextStateNum);
                %                 [NextStateNameArr NextStateNumArr NextStateProbArr] = obj.TransitionSimulateNum(ByteCluster);
                [NextStateNameArr NextStateNumArr NextStateProbArr] = obj.TransitionSimulateNum(NextStateName);
                P = cumsum(NextStateProbArr);
                ii = min(find(P-rand>0));
                if(length(ii) == 0)
                    
                    NextStateNameArr
                end
%                 NextStateNameArr{1}
%                 ii
                NextStateName = NextStateNameArr{ii};
                %                 NextStateNum = NextStateNumArr(ii);
                NextStateNum = NextStateNumArr{ii};
                count = count+1;
            end
            EndStateName = NextStateName;
            EndStateNum = NextStateNum;
            ACstates_AgCap = NAgCapturedNum(obj,EndStateName);
            %             ACstates_AgCap = NAgCapturedNum(obj,EndStateNum+1);
            
            %             ACstates_AgCap = NAgCaptured(obj,EndStateNum+1);
        end
        
        function [NextStateName NextStateNum NextStateProb] = UpdateNextState(obj,ClusterState,bit,BCRj,p,NextStateNum,NextStateProb,NextStateName)
            count = length(NextStateProb)+1;
            NextStateNum(count) = obj.ClusterStateNew(ClusterState,BCRj,bit);
            NextStateProb(count) = p;
            NewSt = ClusterState;
            NewSt(BCRj) = bit;
            %             NextStateName{count} = dec2base(NextStateNum(count),obj.StatesNum);
            NextStateName{count} = NewSt;
            
            %                 NewCState = ClusterState;
            %             NewCState(idx) = NewS;
            %             NewCStateDec = base2dec(NewCState,obj.StatesNum);
            
        end
        
        
        %%%%
        %%%%
        %%%%
        
        function [NextStateName NextStateNum NextStateProb] = UpdateNextStateNum(obj,ClusterState,bit,BCRj,p,NextStateNum,NextStateProb,NextStateName)
            %                                                               UpdateNextStateNum(ByteCluster,4,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
            count = length(NextStateProb)+1;
            NextStateName1 = ClusterState;
            NextStateName1(BCRj) = bit;
            NextStateName{count} = NextStateName1;
            NBCR = obj.NBCR;
            NextStateNum1 = 0;
            for i=1:NBCR
                NextStateNum1 = NextStateNum1+NextStateName1(i)*9^(NBCR-i);
            end
            NextStateNum1 = NextStateNum1+1;
            NextStateNum{count} = NextStateNum1;
            %             NextStateNum(count) = obj.ClusterStateNew(ClusterState,BCRj,bit);
            %             NewCStateDec = ClusterStateNew(obj,ClusterState,idx,NewS);
            %             NewCState = ClusterState;
            %             NewCState(idx) = NewS;
            %             NewCStateDec = base2dec(NewCState,obj.StatesNum);
            
            
            NextStateProb{count} = p;
            %             NewSt = ClusterState;
            %             NewSt(BCRj) = bit;
            %             %             NextStateName{count} = dec2base(NextStateNum(count),obj.StatesNum);
            %             NextStateName{count} = NewSt;
            
            %                 NewCState = ClusterState;
            %             NewCState(idx) = NewS;
            %             NewCStateDec = base2dec(NewCState,obj.StatesNum);
            
        end
        
        
        %%%%
        %%%%
        %%%%
        %
        %         function [NextStateName NextStateNum NextStateProb] = TransitionSimulate(obj,ClusterState,NextStateNum);
        %
        %             NAg = obj.NAg;
        %             kb = 1.38064852e-23;
        %             T=273;
        %             xb = 1e-9;
        %
        %             NBCR = obj.NBCR;
        %             F = obj.F;
        %             % Each BCR is a word (16 bit). The transition states are
        %             Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
        %             Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
        %             TCstates = []; % Transition cluster states
        %             ACstates = []; % Absorbing cluster states
        %             %             StatesNum = length(Tstate)+length(Astate);
        %             %             obj.StatesNum = StatesNum;
        %
        %             %             W = sparse(StatesNum^NBCR,StatesNum^NBCR);
        %             %             W = sparse(1,StatesNum^NBCR);
        %
        %
        %
        %             %             for i=0:(StatesNum^NBCR-1)
        %             %                 ClusterState = dec2base(i,StatesNum,NBCR);
        %
        %             TransitionIndex = regexp(ClusterState,'[0-4]');
        %             TransitionNonConIndex = regexp(ClusterState,'[5-7]');
        %             NumConnected = length(TransitionIndex);
        %             NumNonConnected = length(TransitionNonConIndex);
        %
        %
        %             NBCR = obj.NBCR;
        %             ByteCluster = zeros(1,NBCR);
        %             ClusterState1 = dec2base(NextStateNum,obj.StatesNum,obj.NBCR)
        %             b=NextStateNum;
        %             for p=1:NBCR
        %                 ByteCluster(NBCR-p+1) = mod(b,9);
        %                 b = floor(b/9);
        %             end
        %             TransitionIndex = find(ByteCluster<=4);
        %             TransitionNonConIndex = find( (ByteCluster>=5) & (ByteCluster<7) );
        %             NumConnected1 = length(TransitionIndex);
        %             NumNonConnected1 = length(TransitionNonConIndex);
        %
        %             %             if(NumConnected)
        %             %                 TCstates = [TCstates , i+1];
        %             %             else
        %             %                 ACstates = [ACstates , i+1];
        %             %             end
        %             %Connected states:
        %             f = F/NumConnected;
        %
        %             lambda = exp(-obj.EAgMem)*exp(xb*f/(kb*T));
        %             r = exp(-obj.EAbAg)*exp(xb*f/(kb*T));
        %             k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kb*T));
        %             xi =  exp(-obj.EAgMem)*exp(0.5*xb*f/(kb*T));
        %
        %             nu = 0;
        %
        %             NAgA = obj.NAgAttached(ClusterState);
        %
        %             FreeAg = NAg-NAgA;
        %             if(FreeAg<0)
        %                 FreeAg=0; % cannot be in this state if NAg<NAgA
        %                 %                 continue;
        %             end
        %             q = FreeAg*obj.q;
        %             Z = obj.StatePartition(ClusterState,f,FreeAg,NumConnected);
        %
        %             NextStateNum = [];
        %             NextStateProb = [];
        %             NextStateName = {};
        %
        %             for j=1:NumConnected
        %                 BCRj = TransitionIndex(j);
        %                 BCRState = ClusterState(BCRj);
        %                 switch BCRState
        %                     case '0' % transition to 5, 6, 7
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'5',BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'6',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %                     case '1'
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'5',BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'7',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %
        %                     case '2'
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,nu/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'7',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'8',BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %
        %                     case '3'
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,nu/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'6',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'8',BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %
        %                     case '4'
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'0',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'1',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'2',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
        %                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'3',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %
        %                 end
        %             end
        %             if(NumConnected)
        %                 for j=1:NumNonConnected
        %                     BCRj = TransitionNonConIndex(j);
        %                     BCRState = ClusterState(BCRj);
        %
        %                     switch BCRState
        %                         case '5' % transition to 5, 6, 7
        %                             [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'0',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
        %                             [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'1',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %                         case '6'
        %                             [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'3',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %                         case '7'
        %                             [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'2',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
        %
        %                     end
        %                 end
        %             end
        %             if( round(sum(NextStateProb),2)>1 )
        %                 sum(NextStateProb)
        %             end
        %
        %         end
        %
        %%%%
        %%%%
        %%%%
        function [NextStateName NextStateNum NextStateProb] = TransitionSimulateNum(obj,ClusterState);
            
            NAg = obj.NAg;
            kb = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            NBCR = obj.NBCR;
            F = obj.F;
            % Each BCR is a word (16 bit). The transition states are
            Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            %             StatesNum = length(Tstate)+length(Astate);
            %             obj.StatesNum = StatesNum;
            
            %             W = sparse(StatesNum^NBCR,StatesNum^NBCR);
            %             W = sparse(1,StatesNum^NBCR);
            
            
            
            %             for i=0:(StatesNum^NBCR-1)
            %                 ClusterState = dec2base(i,StatesNum,NBCR);
            
            %             TransitionIndex = regexp(ClusterState,'[0-4]');
            %             TransitionNonConIndex = regexp(ClusterState,'[5-7]');
            %             NumConnected = length(TransitionIndex);
            %             NumNonConnected = length(TransitionNonConIndex);
            
            
            %             NBCR = obj.NBCR;
            %             ByteCluster = zeros(1,NBCR);
            %             ClusterState1 = dec2base(NextStateNum,obj.StatesNum,obj.NBCR)
            %             b=NextStateNum;
            %             for p=1:NBCR
            %                 ByteCluster(NBCR-p+1) = mod(b,9);
            %                 b = floor(i/9);
            %             end
            ByteCluster = ClusterState;
            
            TransitionIndex = find(ByteCluster<=4);
%             TransitionNonConIndex = find( (ByteCluster>=5) & (ByteCluster<7) ); %!!!
            TransitionNonConIndex = find( (ByteCluster>=5) & (ByteCluster<=7) );
            NumConnected = length(TransitionIndex);
            NumNonConnected = length(TransitionNonConIndex);
            
            %             if(NumConnected)
            %                 TCstates = [TCstates , i+1];
            %             else
            %                 ACstates = [ACstates , i+1];
            %             end
            %Connected states:
            f = F/NumConnected;
            
            lambda = exp(-obj.EAgMem)*exp(xb*f/(kb*T));
            r = exp(-obj.EAbAg)*exp(xb*f/(kb*T));
            k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kb*T));
            xi =  exp(-obj.EAgMem)*exp(0.5*xb*f/(kb*T));
            
            nu = 0;
            
            NAgA = obj.NAgAttachedNum(ClusterState);
            
            FreeAg = NAg-NAgA;
            if(FreeAg<0)
                FreeAg=0; % cannot be in this state if NAg<NAgA
                %                 continue;
            end
            q = FreeAg*obj.q;
            Z = obj.StatePartitionNum(ClusterState,f,FreeAg,NumConnected);
            
            NextStateNum = {};
            NextStateProb = {};
            NextStateName = {};
            
            for j=1:NumConnected
                BCRj = TransitionIndex(j);
                %                 BCRState = ClusterState(BCRj);
                BCRState = ByteCluster(BCRj);
                switch BCRState
                    case 0 % transition to 5, 6, 7
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,4,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,5,BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
%                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,6,BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,6,BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);
                        
                    case 1
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,4,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,5,BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
%                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,7,BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,7,BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);                      
                        
                    case 2
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,4,BCRj,nu/Z,NextStateNum,NextStateProb,NextStateName);
%                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,7,BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,7,BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,8,BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);

                    case 3
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,4,BCRj,nu/Z,NextStateNum,NextStateProb,NextStateName);
%                         [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,6,BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,6,BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,8,BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);
                        
                    case 4
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,0,BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,1,BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,2,BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,3,BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        
                end
            end
            if(NumConnected)
                for j=1:NumNonConnected
                    BCRj = TransitionNonConIndex(j);
                    BCRState = ClusterState(BCRj);
                    
                    switch BCRState
                        case 5 % transition to 5, 6, 7
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,0,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,1,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            
                        case 6
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,3,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            
                        case 7
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextStateNum(ByteCluster,2,BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            
                    end
                end
            end
            %             NextStateProb = cell2mat(NextStateProb);
            NextStateProb = [NextStateProb{:}];
            if( round(sum(NextStateProb),2)>1 )
                sum(NextStateProb)
            end
            
        end
        
        %%%
        %%%
        %%%
        %
        %         function NAgA = NAgAttached(obj,ClusterState)
        %             NAgA = 0;
        %             for n=1:obj.NBCR
        %                 BCRState = ClusterState(n);
        %                 switch BCRState
        %                     case '0'
        %                         NAgA = NAgA + 1;
        %                     case '1'
        %                         NAgA = NAgA + 1;
        %                     case '2'
        %                         NAgA = NAgA + 2;
        %                     case '3'
        %                         NAgA = NAgA + 2;
        %                     case '4'
        %                         NAgA = NAgA + 2;
        %                     case '6'
        %                         NAgA = NAgA + 1;
        %                     case '7'
        %                         NAgA = NAgA + 1;
        %                     case '8'
        %                         NAgA = NAgA + 2;
        %                 end
        %             end
        %         end
        %
        %
        %%%
        %%%
        %%%
        
        function NAgA = NAgAttachedNum(obj,ClusterState)
            NAgA = 0;
            for n=1:obj.NBCR
                BCRState = ClusterState(n);
                switch BCRState
                    case 0
                        NAgA = NAgA + 1;
                    case 1
                        NAgA = NAgA + 1;
                    case 2
                        NAgA = NAgA + 2;
                    case 3
                        NAgA = NAgA + 2;
                    case 4
                        NAgA = NAgA + 2;
                    case 6
                        NAgA = NAgA + 1;
                    case 7
                        NAgA = NAgA + 1;
                    case 8
                        NAgA = NAgA + 2;
                end
            end
        end
        
        %%%
        %%%
        %%%
        %
        %         function Z= StatePartition(obj,ClusterState,F,FreeAg,NumConnected)
        %             Z = 0;
        %             NAg = obj.NAg;
        %             kb = 1.38064852e-23;
        %             T=273;
        %             xb = 1e-9;
        %
        %             lambda = exp(-obj.EAgMem)*exp(xb*F/(kb*T));
        %             r = exp(-obj.EAbAg)*exp(xb*F/(kb*T));
        %             k = exp(-obj.EAbAg)*exp(0.5*xb*F/(kb*T));
        %             xi =  exp(-obj.EAgMem)*exp(0.5*xb*F/(kb*T));
        %
        %             NAgA = obj.NAgAttached(ClusterState);
        %             FreeAg = NAg-NAgA;
        %             if(FreeAg<0)
        %                 FreeAg=0; % cannot be in this state if NAg<NAgA
        %             end
        %             q = FreeAg*obj.q;
        %             %             q = (NAg-1)*obj.q;
        %             nu = 0;
        %
        %             for n=1:obj.NBCR
        %                 BCRState = ClusterState(n);
        %                 switch BCRState
        %                     case '0'
        %                         Z = Z + (r+q+xi);
        %                     case '1'
        %                         Z = Z + (r+q+xi);
        %                     case '2'
        %                         Z = Z + (k+nu+lambda);
        %                     case '3'
        %                         Z = Z + (k+nu+lambda);
        %                     case '4'
        %                         Z = Z + (2*k+2*xi);
        %                     case '5'
        %                         if(NumConnected)
        %                             Z = Z + 2*q;
        %                         end
        %                     case '6'
        %                         if(NumConnected)
        %                             Z = Z + q;
        %                         end
        %                     case '7'
        %                         if(NumConnected)
        %                             Z = Z + q;
        %                         end
        %                 end
        %             end
        %         end
        %
        %
        %
        %%%
        %%%
        %%%
        
        
        %%%
        %%%
        %%%
        
        function Z= StatePartitionNum(obj,ClusterState,F,FreeAg,NumConnected)
            Z = 0;
            NAg = obj.NAg;
            kb = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            lambda = exp(-obj.EAgMem)*exp(xb*F/(kb*T));
            r = exp(-obj.EAbAg)*exp(xb*F/(kb*T));
            k = exp(-obj.EAbAg)*exp(0.5*xb*F/(kb*T));
            xi =  exp(-obj.EAgMem)*exp(0.5*xb*F/(kb*T));
            
            NAgA = obj.NAgAttachedNum(ClusterState);
            FreeAg = NAg-NAgA;
            if(FreeAg<0)
                FreeAg=0; % cannot be in this state if NAg<NAgA
            end
            q = FreeAg*obj.q;
            %             q = (NAg-1)*obj.q;
            nu = 0;
            
            for n=1:obj.NBCR
                BCRState = ClusterState(n);
                switch BCRState
                    case 0
%                         Z = Z + (r+q+xi);
                        Z = Z + (r+q+lambda);
                    case 1
%                         Z = Z + (r+q+xi);
                        Z = Z + (r+q+lambda);
                    case 2
%                         Z = Z + (k+nu+lambda);
                        Z = Z + (r+nu+lambda);
                    case 3
%                         Z = Z + (k+nu+lambda);
                        Z = Z + (r+nu+lambda);
                    case 4
                        Z = Z + (2*k+2*xi);
                    case 5
                        if(NumConnected)
                            Z = Z + 2*q;
                        end
                    case 6
                        if(NumConnected)
                            Z = Z + q;
                        end
                    case 7
                        if(NumConnected)
                            Z = Z + q;
                        end
                end
            end
        end
        
        %%%
        %%%
        %%%
        
        
        %         function NewCStateDec = ClusterStateNew(obj,ClusterState,idx,NewS);
        %             NewCState = ClusterState;
        %             NewCState(idx) = NewS;
        %             NewCStateDec = base2dec(NewCState,obj.StatesNum);
        %         end
        
        %%%
        %%%
        %%%
        
        %         function ACstates_AgCap = NAgCaptured(obj,ACstates)
        %
        %             ACstates_AgCap = zeros(1,length(ACstates));
        %             for i=1:length(ACstates_AgCap)
        %                 NAgC = 0;
        %                 ClusterState = dec2base(ACstates(i)-1,obj.StatesNum,obj.NBCR);
        %                 for n=1:size(ClusterState,2)
        %                     %                 ACstates,StatesNum,obj.NBCR)
        %                     BCRState = ClusterState(n);
        %                     switch BCRState
        %                         case '5'
        %                             NAgC = NAgC + 0;
        %                         case '6'
        %                             NAgC = NAgC + 1;
        %                         case '7'
        %                             NAgC = NAgC + 1;
        %                         case '8'
        %                             NAgC = NAgC + 2;
        %                     end
        %                 end
        %                 ACstates_AgCap(i) = NAgC;
        %             end
        %         end
        %
        %%%
        %%%
        %%%
        
        
        %%%
        %%%
        %%%
        
        %         function ACstates_AgCap = NAgCapturedNum(obj,ACstates)
        %         function ACstates_AgCap = NAgCapturedNum(obj,ACstates)
        function NAgC = NAgCapturedNum(obj,ACstate)
            
            %             ACstates_AgCap = zeros(1,size(ACstates,1));
            %             for i=1:length(ACstates_AgCap)
            NAgC = 0;
            ClusterState = ACstate;
            %                 ClusterState = dec2base(ACstates(i)-1,obj.StatesNum,obj.NBCR);
            
            %                 NBCR = obj.NBCR;
            %                 ByteCluster = zeros(1,NBCR);
            %                 b=ACstates-1;
            %                 for p=1:NBCR
            %                     ByteCluster(NBCR-p+1) = mod(b,9);
            %                     b = floor(b/9);
            %                 end
            %                 ClusterState = ByteCluster;
            %                 dd = ByteCluster;
            
            for n=1:size(ClusterState,2)
                %                 ACstates,StatesNum,obj.NBCR)
                BCRState = ClusterState(n);
                switch BCRState
                    case 5
                        NAgC = NAgC + 0;
                    case 6
                        NAgC = NAgC + 1;
                    case 7
                        NAgC = NAgC + 1;
                    case 8
                        NAgC = NAgC + 2;
                end
            end
            %                 ACstates_AgCap(i) = NAgC;
            %             end
        end
        
        %%%
        %%%
        %%%
        
        function Flag = isAbsorbing(obj,State)
            idx = find(State <5 );
            if(length(idx))
                Flag = 0;
            else
                Flag = 1;
            end
        end
        
    end
end