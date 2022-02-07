classdef AbAgCluster < handle
    properties
        NBCR;
        NAg;
        W;
        WnBCR;
        Weq;
        q;
        tau;
        F;
        EAgMem;
        EAbAg;
        TCstates;
        ACstates;
        StatesNum;
        AgCapPro;
        AgCapProArr;
        AgCapProStartingState;
        AgCapProStartingStateArr;
        StartingState;
        StartingStateNames;
        StartingStateEq;
        StartingStateNamesEq;
        Parr;
        ParrNBCR;
        ClusterSize;
        AgCapState;
        TstateName;
        AstateName;
        ICnum;
        aff;
        AgCapProbSum;
        AgCapProbSumFlag;
        SumCapProbMat;
    end
    methods
        function obj = AbAgCluster(aff,NBCR,NAg,q,EAbAg,tau,F,ClusterSize,EqFlag)
            if(nargin > 0)
                
                %                 EqFlag = 1;
                obj.aff = aff;
                obj.WnBCR = {};
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
                else
                end
                obj.EAgMem = aff.EAgMem;
                obj.tau=tau;
                obj.F=F;
                obj.AgCapProbSumFlag = 0;
                switch aff.BindFlag
                    case 8
                        obj.ConstructWMat();
                        obj.CaptureProb();
                        obj.StartingStatesFind(obj.TCstates);
                    case 9
                        obj.ConstructWMat2Ag();
                        obj.CaptureProb2Ag();
                        obj.StartingStatesFind2Ag(obj.TCstates);
                    case 10
                        obj.ConstructWMat2AgExplicitKon();
                        obj.CaptureProb2Ag();
                        obj.StartingStatesFind2Ag(obj.TCstates);
                end
                if(EqFlag)
                    obj.ConstructEQWMat();
                end
                
                if(EqFlag==0)
                    switch aff.BindFlag
                        case 8
                            obj.CalcProbDensityNoneq();
                            [C,ia,ib] = intersect(obj.TCstates,obj.StartingState);
                            %                     obj.AgCapProStartingState = obj.AgCapPro(obj.StartingState,:);
                            obj.AgCapProStartingState = obj.AgCapPro(ia,:);
                        case 9
                            obj.CalcProbDensityNoneq2Ag();
                            [C,ia,ib] = intersect(obj.TCstates,obj.StartingState);
%                             obj.AgCapProStartingState = obj.AgCapPro(ia,:);
                            obj.AgCapProStartingStateArr = cellfun(@(x)x(ia,:),obj.AgCapProArr,'un',0);
                        case 10
                            obj.CalcProbDensityNoneq2AgExplicitDensity();
                            [C,ia,ib] = intersect(obj.TCstates,obj.StartingState);
%                             obj.AgCapProStartingState = obj.AgCapPro(ia,:);
                            obj.AgCapProStartingStateArr = cellfun(@(x)x(ia,:),obj.AgCapProArr,'un',0);
                            
                    end
                else
                    [C,ia,ib] = intersect(obj.TCstates,obj.StartingStateEq);
                    obj.AgCapProStartingState = obj.AgCapPro(ia,:);
                end
                
                %                 obj.AbsorbingStatesFind();
                %                 obj.StartingStatesFind(obj.TCstates);
                %                 obj.CalcProbDensityNoneq();
                
            else
            end
        end
        
        
        %%%%%
        %%%%%
        %%%%%
        
        function ConstructWMat2Ag(obj)
            
            global ICArr;
            %             NState = 2^(2*obj.NBCR);
            NAg = obj.NAg;
            kB = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            NBCR = obj.NBCR;
            %             F = 1e-12;
            F = obj.F;
            %             F=0;
            
            
            % Each BCR is a word (16 bit). The transition states are
            Tstate = [121,202,283,337,355,364,418,444, 445,526, 598, 604, 607, 685, 687, 688];
            TstateName = { '011111', '021111', '101111', '110111', '111011', '111111', '120111', '121110', '121111', '201111', '211011', '211101', '211111', '221101', '221110', '221111'};  % connected states : force falls on the bond
            obj.TstateName = TstateName;
            Astate = [40,112,201,256,328,417,523,595,684];
            AstateName = {'001111', '011011', '021110', '100111', '110011', '120110', '201101', '211001', '221100'};
            obj.AstateName = AstateName;
            statesNames = [TstateName AstateName];
            % 001111 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            
            % Each BCR is a word (16 bit). The transition states are
            %             Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            %             Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            StatesNum = length(Tstate)+length(Astate);
            obj.StatesNum = StatesNum;
            
            TotICNum = length(ICArr{1})+length(ICArr{2});
            if( rand< length(ICArr{1})/TotICNum)
                itype = 1;
            else
                itype = 2;
            end
            iIC = randi([1 length(ICArr{itype})]);
%             obj.ICnum = [1 1];
            obj.ICnum = [itype iIC];
            IC = ICArr{obj.ICnum(1)}{obj.ICnum(2)};
            %             NAga = IC(1);
            %             NAgb = IC(2);
            AgAttachMat = zeros(StatesNum^NBCR,2);
            for i=0:(StatesNum^NBCR-1)
                [NAgaA NAgbA] = obj.NAgAttached2AgInt(i);
                AgAttachMat(i+1,1) = NAgaA;
                AgAttachMat(i+1,2) = NAgbA;
            end
            
            for n_bcr=1:obj.aff.NBCR
                
                NAga = IC{n_bcr}(1);
                NAgb = IC{n_bcr}(2);
                W = sparse(StatesNum^NBCR,StatesNum^NBCR);
                for i=0:(StatesNum^NBCR-1)
%                     ClusterState = num2str(i);
              
                    if(i<16)
                        NumConnected = 1;
                        NumNonConnected = 0;
                        if(n_bcr==1)
                            TCstates = [TCstates , i+1];
                        end
                    else
                        NumConnected = 1;
                        NumNonConnected = 0;
                        if(n_bcr==1)
                            ACstates = [ACstates , i+1];
                        end
                    end

                    %                 if(NumConnected)
                    %                     TCstates = [TCstates , i+1];
                    %                 else
                    %                     ACstates = [ACstates , i+1];
                    %                 end
                    %Connected states:
                    f = F/NumConnected;
                    %                 NumConnected
                    lambda = exp(-obj.EAgMem)*exp(xb*f/(kB*T));
                    r = exp(-obj.EAbAg)*exp(xb*f/(kB*T));
                    k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kB*T));
                    xi =  exp(-obj.EAgMem)*exp(0.5*xb*f/(kB*T));
                    
                    
%                     EAbAg_a = obj.aff.wcv(1)*obj.EAbAg(1)+obj.aff.wcv(2)*obj.EAbAg(2);
%                     EAbAg_b = obj.aff.wcv(1)*obj.EAbAg(1)+obj.aff.wcv(2)*obj.EAbAg(3);
%                     EAbAg_a = obj.aff.rcv(1)*obj.EAbAg(1)+obj.aff.rcv(2)*obj.EAbAg(2);
%                     EAbAg_b = obj.aff.rcv(1)*obj.EAbAg(1)+obj.aff.rcv(2)*obj.EAbAg(3);
%                     ea = obj.EAbAg(2)/(obj.EAbAg(2)+obj.EAbAg(3));
%                     eb = obj.EAbAg(3)/(obj.EAbAg(2)+obj.EAbAg(3));
%                     if(isnan(ea) | isnan(eb))
%                         ea=0;
%                         eb=0;
%                     end
%                     ea=1;eb=1;
%                     EAbAg_a = obj.aff.wcv(1)*obj.EAbAg(1)+obj.aff.wcv(2)*ea*obj.EAbAg(2);
%                     EAbAg_b = obj.aff.wcv(1)*obj.EAbAg(1)+obj.aff.wcv(2)*eb*obj.EAbAg(3);

                    EAbAg_a = obj.EAbAg(1)+obj.EAbAg(2);
                    EAbAg_b = obj.EAbAg(1)+obj.EAbAg(3);
%                     r = exp(-obj.EAbAg)*exp(xb*f/(kB*T));
%                     k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kB*T));
                    
%                     ra=1;
%                     rb=1;
%                     ka=1;
%                     kb=1;
%                     if(EAbAg_a<0)
%                         EAbAg_a
%                     end
                    ra = exp(-EAbAg_a)*exp(xb*f/(kB*T));
                    rb = exp(-EAbAg_b)*exp(xb*f/(kB*T));
                    ka = exp(-EAbAg_a)*exp(0.5*xb*f/(kB*T));
                    kb = exp(-EAbAg_b)*exp(0.5*xb*f/(kB*T));
                    
                    nu = 0;
                    
%                     [NAgaA NAgbA] = obj.NAgAttached2Ag(ClusterState);
%                     [NAgaA NAgbA] = obj.NAgAttached2AgInt(i);
%                     FreeAga = NAga-NAgaA;
%                     FreeAgb = NAgb-NAgbA;
                    
                    FreeAga = NAga-AgAttachMat(i+1,1);
                    FreeAgb = NAgb-AgAttachMat(i+1,2);
                    
                    %                 if(FreeAg<0)
                    if( (FreeAga<0) | (FreeAgb<0) )
                        %                     FreeAg=0; % cannot be in this state if NAg<NAgA
                        continue;
                    end
                    %                 q = FreeAg*obj.q;
                    
%                     qa0 = obj.aff.wcv(1)*obj.q(1)+obj.aff.wcv(2)*obj.q(2);
%                     qb0 = obj.aff.wcv(1)*obj.q(1)+obj.aff.wcv(2)*obj.q(3);
%                     qa0 = obj.aff.rcv(1)*obj.q(1)+obj.aff.rcv(2)*obj.q(2);
%                     qb0 = obj.aff.rcv(1)*obj.q(1)+obj.aff.rcv(2)*obj.q(3);
                    
%                     ea = obj.q(2)/(obj.q(2)+obj.q(3));
%                     eb = obj.q(3)/(obj.q(2)+obj.q(3));
%                     if(isnan(ea) | isnan(eb))
%                         ea=0;
%                         eb=0;
%                     end
%                     ea=1;eb=1;
%                     qa0 = obj.aff.wcv(1)*obj.q(1)+obj.aff.wcv(2)*ea*obj.q(2);
%                     qb0 = obj.aff.wcv(1)*obj.q(1)+obj.aff.wcv(2)*eb*obj.q(3);
                    
                    qa0 = obj.q(1)+obj.q(2);
                    qb0 = obj.q(1)+obj.q(3);
%                     qa=FreeAga*obj.q;
%                     qb=FreeAgb*obj.q;
                    qa=FreeAga*qa0;
                    qb=FreeAgb*qb0;
                    
                    %                 q = (NAg-1)*obj.q;
                    %                 Z = obj.StatePartition(ClusterState,f,FreeAg,NumConnected);
                    for j=1:NumConnected
                        %                     BCRj = TransitionIndex(j);
                        %                     BCRState = ClusterState(BCRj);
                        %                     BCRState = statesNames{i+1};
%                         BCRState = num2str(i);
                        BCRState = i;
                        switch BCRState
                            case 0 % transition to 5, 6, 7
                                Z = qa+qb+ra+lambda;
                                
                                p0_5 = qa/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p0_5;
                                
                                p0_12 = qb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'12');
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = p0_12;
                                
                                p0_16 = ra/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = p0_16;
                                
                                %                             p06 = xi/Z;
                                p0_17 = lambda/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'17');
                                NewCStateDec = 17;
                                W(i+1,NewCStateDec+1) = p0_17;
                            case 1
                                Z = qa+qb+rb+lambda;
                                
                                %                             p_8 = qa/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = qa/Z;
                                
                                %                             p15 = qb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'15');
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = qb/Z;
                                
                                %                             p16 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = rb/Z;
                                
                                %                             p18 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'18');
                                NewCStateDec = 18;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 2
                                Z = qa+qb+ra+lambda;
                                
                                p2_5 = qa/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p2_5;
                                
                                p2_8 = qb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = p2_8;
                                
                                p2_16 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = p2_16;
                                
                                p2_19 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'19');
                                NewCStateDec = 19;
                                W(i+1,NewCStateDec+1) = p2_19;
                                
                                %                             p24 = nu/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                                %                             W(i+1,NewCStateDec+1) = p24;
                                %
                                % %                             p27 = k/Z;
                                %                             p27 = r/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                                %                             W(i+1,NewCStateDec+1) = p27;
                                %
                                %                             p28 = lambda/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                %                             W(i+1,NewCStateDec+1) = p28;
                                
                            case 3
                                Z = nu+ra+lambda;
                                
                                p3_5 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p3_5;
                                
                                p3_19 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'19');
                                NewCStateDec = 19;
                                W(i+1,NewCStateDec+1) = p3_19;
                                
                                p3_20 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'20');
                                NewCStateDec = 20;
                                W(i+1,NewCStateDec+1) = p3_20;
                                
                            case 4
                                Z = nu+ra+lambda;
                                
                                p3_5 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p3_5;
                                
                                p3_17 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'17');
                                NewCStateDec = 17;
                                W(i+1,NewCStateDec+1) = p3_17;
                                
                                p3_20 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'20');
                                NewCStateDec = 20;
                                W(i+1,NewCStateDec+1) = p3_20;
                                
                            case 5
                                Z = 2*ka+2*xi;
                                
                                p5_0 = ka/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                                NewCStateDec = 0;
                                W(i+1,NewCStateDec+1) = p5_0;
                                
                                p5_2 = ka/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                                NewCStateDec = 2;
                                W(i+1,NewCStateDec+1) = p5_2;
                                
                                p5_3 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                                NewCStateDec = 3;
                                W(i+1,NewCStateDec+1) = p5_3;
                                
                                p5_4 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                                NewCStateDec = 4;
                                W(i+1,NewCStateDec+1) = p5_4;
                                
                            case 6
                                Z = nu+rb+lambda;
                                
                                p6_8 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = p6_8;
                                
                                p6_19 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'19');
                                NewCStateDec = 19;
                                W(i+1,NewCStateDec+1) = p6_19;
                                
                                p6_21 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'21');
                                NewCStateDec = 21;
                                W(i+1,NewCStateDec+1) = p6_21;
                                
                            case 7
                                Z = nu+ra+lambda;
                                
                                p7_8 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = p7_8;
                                
                                p7_18 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'18');
                                NewCStateDec = 18;
                                W(i+1,NewCStateDec+1) = p7_18;
                                
                                p7_21 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'21');
                                NewCStateDec = 21;
                                W(i+1,NewCStateDec+1) = p7_21;
                                
                            case 8
                                Z = ka+kb+2*xi;
                                
                                p8_1 = ka/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                                NewCStateDec = 1;
                                W(i+1,NewCStateDec+1) = p8_1;
                                
                                p8_2 = kb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                                NewCStateDec = 2;
                                W(i+1,NewCStateDec+1) = p8_2;
                                
                                p8_6 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                                NewCStateDec = 6;
                                W(i+1,NewCStateDec+1) = p8_6;
                                
                                p8_7 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                                NewCStateDec = 7;
                                W(i+1,NewCStateDec+1) = p8_7;
                                
                            case 9
                                Z = qa+qb+rb+lambda;
                                
                                p9_12 = qa/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'12');
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = p9_12;
                                
                                p9_15 = qb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'15');
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = p9_15;
                                
                                p9_16 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = p9_16;
                                
                                p9_22 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'22');
                                NewCStateDec = 22;
                                W(i+1,NewCStateDec+1) = p9_22;
                                
                            case 10
                                Z = nu+rb+lambda;
                                
                                p10_12 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'12');
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = p10_12;
                                
                                p7_17 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'17');
                                NewCStateDec = 17;
                                W(i+1,NewCStateDec+1) = p7_17;
                                
                                p7_23 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'23');
                                NewCStateDec = 23;
                                W(i+1,NewCStateDec+1) = p7_23;
                                
                            case 11
                                Z = nu+ra+lambda;
                                
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = nu/Z;
                                
                                NewCStateDec = 22;
                                W(i+1,NewCStateDec+1) = ra/Z;
                                
                                NewCStateDec = 23;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 12
                                Z = kb+ka+2*xi;
                                
                                NewCStateDec = 0;
                                W(i+1,NewCStateDec+1) = kb/Z;
                                
                                NewCStateDec = 9;
                                W(i+1,NewCStateDec+1) = ka/Z;
                                
                                NewCStateDec = 10;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                                NewCStateDec = 11;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                            case 13
                                Z = nu+rb+lambda;
                                
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = nu/Z;
                                
                                NewCStateDec = 22;
                                W(i+1,NewCStateDec+1) = rb/Z;
                                
                                NewCStateDec = 24;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 14
                                Z = nu+rb+lambda;
                                
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = nu/Z;
                                
                                NewCStateDec = 18;
                                W(i+1,NewCStateDec+1) = rb/Z;
                                
                                NewCStateDec = 24;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 15
                                Z = 2*kb+2*xi;
                                
                                NewCStateDec = 1;
                                W(i+1,NewCStateDec+1) = kb/Z;
                                
                                NewCStateDec = 9;
                                W(i+1,NewCStateDec+1) = kb/Z;
                                
                                NewCStateDec = 13;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                                NewCStateDec = 14;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                        end
                    end
                    %                 if(NumConnected)
                    %                     for j=1:NumNonConnected
                    %                         BCRj = TransitionNonConIndex(j);
                    %                         BCRState = ClusterState(BCRj);
                    %
                    % %                         switch BCRState
                    % %                             case '5' % transition to 5, 6, 7
                    % %                                 p50 = q/Z;
                    % %                                 % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                    % %                                 W(i+1,NewCStateDec+1) = p50;
                    % %
                    % %                                 p51 = q/Z;
                    % %                                 % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                    % %                                 W(i+1,NewCStateDec+1) = p51;
                    % %                             case '6'
                    % %                                 p63 = q/Z;
                    % %                                 % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                    % %                                 W(i+1,NewCStateDec+1) = p63;
                    % %                             case '7'
                    % %                                 p72 = q/Z;
                    % %                                 % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                    % %                                 W(i+1,NewCStateDec+1) = p72;
                    % %                         end
                    %                     end
                    %                 end
                    
%                     if( round(sum(W(i+1,:)),2)>1 )
%                         sum(W(i+1,:));
%                     end
                end
                obj.W = W;
                obj.WnBCR{n_bcr} = W;
            end
            obj.TCstates = TCstates;
            obj.ACstates = ACstates;
            
        end
        
        %%%%%
        %%%%%
        %%%%%
      
        
        %%%%%
        %%%%%  Antigen capture probabilities when we have the explicit
        %%%%%  dependence on the Ag density
        %%%%%
        
        
        function ConstructWMat2AgExplicitKon(obj)
            
            global ICArr;
            %             NState = 2^(2*obj.NBCR);
            NAg = obj.NAg;
            kB = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            NBCR = obj.NBCR;
            %             F = 1e-12;
            F = obj.F;
            %             F=0;
            
            
            % Each BCR is a word (16 bit). The transition states are
            Tstate = [121,202,283,337,355,364,418,444, 445,526, 598, 604, 607, 685, 687, 688];
            TstateName = { '011111', '021111', '101111', '110111', '111011', '111111', '120111', '121110', '121111', '201111', '211011', '211101', '211111', '221101', '221110', '221111'};  % connected states : force falls on the bond
            obj.TstateName = TstateName;
            Astate = [40,112,201,256,328,417,523,595,684];
            AstateName = {'001111', '011011', '021110', '100111', '110011', '120110', '201101', '211001', '221100'};
            obj.AstateName = AstateName;
            statesNames = [TstateName AstateName];
            % 001111 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            
            % Each BCR is a word (16 bit). The transition states are
            %             Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            %             Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            StatesNum = length(Tstate)+length(Astate);
            obj.StatesNum = StatesNum;
            
            TotICNum = length(ICArr{1})+length(ICArr{2});
            if( rand< length(ICArr{1})/TotICNum)
                itype = 1;
            else
                itype = 2;
            end
            iIC = randi([1 length(ICArr{itype})]);
%             obj.ICnum = [1 1];
            obj.ICnum = [itype iIC];
            IC = ICArr{obj.ICnum(1)}{obj.ICnum(2)};
            %             NAga = IC(1);
            %             NAgb = IC(2);
            AgAttachMat = zeros(StatesNum^NBCR,2);
            for i=0:(StatesNum^NBCR-1)
                [NAgaA NAgbA] = obj.NAgAttached2AgInt(i);
                AgAttachMat(i+1,1) = NAgaA;
                AgAttachMat(i+1,2) = NAgbA;
            end
            
            for n_bcr=1:obj.aff.NBCR
                
                NAga = IC{n_bcr}(1);
                NAgb = IC{n_bcr}(2);
                W = sparse(StatesNum^NBCR,StatesNum^NBCR);
                for i=0:(StatesNum^NBCR-1)
%                     ClusterState = num2str(i);
              
                    if(i<16)
                        NumConnected = 1;
                        NumNonConnected = 0;
                        if(n_bcr==1)
                            TCstates = [TCstates , i+1];
                        end
                    else
                        NumConnected = 1;
                        NumNonConnected = 0;
                        if(n_bcr==1)
                            ACstates = [ACstates , i+1];
                        end
                    end

                    f = F/NumConnected;

                    lambda = exp(-obj.EAgMem)*exp(xb*f/(kB*T));
                    r = exp(-obj.EAbAg)*exp(xb*f/(kB*T));
                    k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kB*T));
                    xi =  exp(-obj.EAgMem)*exp(0.5*xb*f/(kB*T));
     
                    EAbAg_a = obj.EAbAg(1)+obj.EAbAg(2);
                    EAbAg_b = obj.EAbAg(1)+obj.EAbAg(3);

                    ra = exp(-EAbAg_a)*exp(xb*f/(kB*T));
                    rb = exp(-EAbAg_b)*exp(xb*f/(kB*T));
                    ka = exp(-EAbAg_a)*exp(0.5*xb*f/(kB*T));
                    kb = exp(-EAbAg_b)*exp(0.5*xb*f/(kB*T));
                    
                    nu = 0;
        
                    FreeAga = NAga-AgAttachMat(i+1,1);
                    FreeAgb = NAgb-AgAttachMat(i+1,2);
                    
                    %                 if(FreeAg<0)
                    if( (FreeAga<0) | (FreeAgb<0) )
                        %                     FreeAg=0; % cannot be in this state if NAg<NAgA
                        continue;
                    end
  
                    qa0 = obj.q(1)+obj.q(2);
                    qb0 = obj.q(1)+obj.q(3);
%                     qa=FreeAga*obj.q;
%                     qb=FreeAgb*obj.q;
%                     Epitope = 1;
%                     qa1Arm = DensityDependentRate(NAga,obj.aff.Epitope,1); % binding rate to Ag "a" for the first arm
%                     qb1Arm = DensityDependentRate(NAgb,obj.aff.Epitope,1); % binding rate to Ag "b" for the first arm
%                     qa2ndArm = DensityDependentRate(NAga,obj.aff.Epitope,2); % binding rate to Ag "a" for the second arm given that the first one is bound
%                     qb2ndArm = DensityDependentRate(NAgb,obj.aff.Epitope,2); % binding rate to Ag "b" for the second arm given that the first one is bound


                    qa1Arm = obj.DensityDependentRate1Arm(NAga,obj.aff.Epitope); % binding rate to Ag "a" for the first arm
                    qb1Arm = obj.DensityDependentRate1Arm(NAgb,obj.aff.Epitope); % binding rate to Ag "b" for the first arm
                    qa2ndArm = obj.DensityDependentRate2Arm(NAga,obj.aff.Epitope); % binding rate to Ag "a" for the second arm given that the first one is bound
                    qb2ndArm = obj.DensityDependentRate2Arm(NAgb,obj.aff.Epitope); % binding rate to Ag "b" for the second arm given that the first one is bound

        

                    
%                     qa=FreeAga*qa0;
%                     qb=FreeAgb*qb0;
                    
                    for j=1:NumConnected

                        BCRState = i;
                        switch BCRState
                            case 0 % transition to 5, 6, 7
%                                 Z = qa+qb+ra+lambda;
                                Z = qa1Arm+qb1Arm+ra+lambda;
                                
%                                 p0_5 = qa/Z;
                                p0_5 = qa1Arm/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p0_5;
                                
%                                 p0_12 = qb/Z;
                                p0_12 = qb1Arm/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'12');
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = p0_12;
                                
                                p0_16 = ra/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = p0_16;
                                
                                %                             p06 = xi/Z;
                                p0_17 = lambda/Z;
                                %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'17');
                                NewCStateDec = 17;
                                W(i+1,NewCStateDec+1) = p0_17;
                            case 1
%                                 Z = qa+qb+rb+lambda;  % Kon2 (on rate of the second arm)
                                Z = qa2ndArm+qb2ndArm+rb+lambda;  % Kon2 (on rate of the second arm)
                                
                                %                             p_8 = qa/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
%                                 W(i+1,NewCStateDec+1) = qa/Z;
                                W(i+1,NewCStateDec+1) = qa2ndArm/Z;
                                
                                %                             p15 = qb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'15');
                                NewCStateDec = 15;
%                                 W(i+1,NewCStateDec+1) = qb/Z;
                                W(i+1,NewCStateDec+1) = qb2ndArm/Z;
                                
                                %                             p16 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = rb/Z;
                                
                                %                             p18 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'18');
                                NewCStateDec = 18;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 2
%                                 Z = qa+qb+ra+lambda; % Kon2 (on rate of the second are)
                                Z = qa2ndArm+qb2ndArm+ra+lambda; % Kon2 (on rate of the second are)
                                
%                                 p2_5 = qa/Z;
                                p2_5 = qa2ndArm/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p2_5;
                                
%                                 p2_8 = qb/Z;
                                p2_8 = qb2ndArm/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = p2_8;
                                
                                p2_16 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = p2_16;
                                
                                p2_19 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'19');
                                NewCStateDec = 19;
                                W(i+1,NewCStateDec+1) = p2_19;
                                                                
                            case 3
                                Z = nu+ra+lambda;
                                
                                p3_5 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p3_5;
                                
                                p3_19 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'19');
                                NewCStateDec = 19;
                                W(i+1,NewCStateDec+1) = p3_19;
                                
                                p3_20 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'20');
                                NewCStateDec = 20;
                                W(i+1,NewCStateDec+1) = p3_20;
                                
                            case 4
                                Z = nu+ra+lambda;
                                
                                p3_5 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                                NewCStateDec = 5;
                                W(i+1,NewCStateDec+1) = p3_5;
                                
                                p3_17 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'17');
                                NewCStateDec = 17;
                                W(i+1,NewCStateDec+1) = p3_17;
                                
                                p3_20 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'20');
                                NewCStateDec = 20;
                                W(i+1,NewCStateDec+1) = p3_20;
                                
                            case 5
                                Z = 2*ka+2*xi;
                                
                                p5_0 = ka/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                                NewCStateDec = 0;
                                W(i+1,NewCStateDec+1) = p5_0;
                                
                                p5_2 = ka/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                                NewCStateDec = 2;
                                W(i+1,NewCStateDec+1) = p5_2;
                                
                                p5_3 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                                NewCStateDec = 3;
                                W(i+1,NewCStateDec+1) = p5_3;
                                
                                p5_4 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                                NewCStateDec = 4;
                                W(i+1,NewCStateDec+1) = p5_4;
                                
                            case 6
                                Z = nu+rb+lambda;
                                
                                p6_8 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = p6_8;
                                
                                p6_19 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'19');
                                NewCStateDec = 19;
                                W(i+1,NewCStateDec+1) = p6_19;
                                
                                p6_21 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'21');
                                NewCStateDec = 21;
                                W(i+1,NewCStateDec+1) = p6_21;
                                
                            case 7
                                Z = nu+ra+lambda;
                                
                                p7_8 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                                NewCStateDec = 8;
                                W(i+1,NewCStateDec+1) = p7_8;
                                
                                p7_18 = ra/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'18');
                                NewCStateDec = 18;
                                W(i+1,NewCStateDec+1) = p7_18;
                                
                                p7_21 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'21');
                                NewCStateDec = 21;
                                W(i+1,NewCStateDec+1) = p7_21;
                                
                            case 8
                                Z = ka+kb+2*xi;
                                
                                p8_1 = ka/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                                NewCStateDec = 1;
                                W(i+1,NewCStateDec+1) = p8_1;
                                
                                p8_2 = kb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                                NewCStateDec = 2;
                                W(i+1,NewCStateDec+1) = p8_2;
                                
                                p8_6 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                                NewCStateDec = 6;
                                W(i+1,NewCStateDec+1) = p8_6;
                                
                                p8_7 = xi/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                                NewCStateDec = 7;
                                W(i+1,NewCStateDec+1) = p8_7;
                                
                            case 9
%                                 Z = qa+qb+rb+lambda;
                                Z = qa2ndArm+qb2ndArm+rb+lambda;
                                
                                
%                                 p9_12 = qa/Z;
                                p9_12 = qa2ndArm/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'12');
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = p9_12;
                                
%                                 p9_15 = qb/Z;
                                p9_15 = qb2ndArm/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'15');
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = p9_15;
                                
                                p9_16 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'16');
                                NewCStateDec = 16;
                                W(i+1,NewCStateDec+1) = p9_16;
                                
                                p9_22 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'22');
                                NewCStateDec = 22;
                                W(i+1,NewCStateDec+1) = p9_22;
                                
                            case 10
                                Z = nu+rb+lambda;
                                
                                p10_12 = nu/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'12');
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = p10_12;
                                
                                p7_17 = rb/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'17');
                                NewCStateDec = 17;
                                W(i+1,NewCStateDec+1) = p7_17;
                                
                                p7_23 = lambda/Z;
                                % NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'23');
                                NewCStateDec = 23;
                                W(i+1,NewCStateDec+1) = p7_23;
                                
                            case 11
                                Z = nu+ra+lambda;
                                
                                NewCStateDec = 12;
                                W(i+1,NewCStateDec+1) = nu/Z;
                                
                                NewCStateDec = 22;
                                W(i+1,NewCStateDec+1) = ra/Z;
                                
                                NewCStateDec = 23;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 12
                                Z = kb+ka+2*xi;
                                
                                NewCStateDec = 0;
                                W(i+1,NewCStateDec+1) = kb/Z;
                                
                                NewCStateDec = 9;
                                W(i+1,NewCStateDec+1) = ka/Z;
                                
                                NewCStateDec = 10;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                                NewCStateDec = 11;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                            case 13
                                Z = nu+rb+lambda;
                                
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = nu/Z;
                                
                                NewCStateDec = 22;
                                W(i+1,NewCStateDec+1) = rb/Z;
                                
                                NewCStateDec = 24;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 14
                                Z = nu+rb+lambda;
                                
                                NewCStateDec = 15;
                                W(i+1,NewCStateDec+1) = nu/Z;
                                
                                NewCStateDec = 18;
                                W(i+1,NewCStateDec+1) = rb/Z;
                                
                                NewCStateDec = 24;
                                W(i+1,NewCStateDec+1) = lambda/Z;
                                
                            case 15
                                Z = 2*kb+2*xi;
                                
                                NewCStateDec = 1;
                                W(i+1,NewCStateDec+1) = kb/Z;
                                
                                NewCStateDec = 9;
                                W(i+1,NewCStateDec+1) = kb/Z;
                                
                                NewCStateDec = 13;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                                NewCStateDec = 14;
                                W(i+1,NewCStateDec+1) = xi/Z;
                                
                        end
                    end

                end
                obj.W = W;
                obj.WnBCR{n_bcr} = W;
            end
            obj.TCstates = TCstates;
            obj.ACstates = ACstates;
            
        end
        
 
        
        
        %%%%%
        %%%%%
        %%%%%
        
        
        
        function ConstructWMat(obj)
            %             NState = 2^(2*obj.NBCR);
            NAg = obj.NAg;
            kb = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            NBCR = obj.NBCR;
            %             F = 1e-12;
            F = obj.F;
            %             F=0;
            % Each BCR is a word (16 bit). The transition states are
            Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            StatesNum = length(Tstate)+length(Astate);
            obj.StatesNum = StatesNum;
            W = sparse(StatesNum^NBCR,StatesNum^NBCR);
            for i=0:(StatesNum^NBCR-1)
                ClusterState = dec2base(i,StatesNum,NBCR);
                
                TransitionIndex = regexp(ClusterState,'[0-4]');
                TransitionNonConIndex = regexp(ClusterState,'[5-7]');
                NumConnected = length(TransitionIndex);
                NumNonConnected = length(TransitionNonConIndex);
                if(NumConnected)
                    TCstates = [TCstates , i+1];
                else
                    ACstates = [ACstates , i+1];
                end
                %Connected states:
                f = F/NumConnected;
                
                lambda = exp(-obj.EAgMem)*exp(xb*f/(kb*T));
                r = exp(-obj.EAbAg)*exp(xb*f/(kb*T));
                k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kb*T));
                xi =  exp(-obj.EAgMem)*exp(0.5*xb*f/(kb*T));
                
                nu = 0;
                
                NAgA = obj.NAgAttached(ClusterState);
                
                FreeAg = NAg-NAgA;
                if(FreeAg<0)
                    FreeAg=0; % cannot be in this state if NAg<NAgA
                    continue;
                end
                q = FreeAg*obj.q;
                %                 q = (NAg-1)*obj.q;
                Z = obj.StatePartition(ClusterState,f,FreeAg,NumConnected);
                for j=1:NumConnected
                    BCRj = TransitionIndex(j);
                    BCRState = ClusterState(BCRj);
                    switch BCRState
                        case '0' % transition to 5, 6, 7
                            p04 = q/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            W(i+1,NewCStateDec+1) = p04;
                            
                            p05 = r/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                            W(i+1,NewCStateDec+1) = p05;
                            
                            %                             p06 = xi/Z;
                            p06 = lambda/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                            W(i+1,NewCStateDec+1) = p06;
                        case '1'
                            p14 = q/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            W(i+1,NewCStateDec+1) = p14;
                            
                            p15 = r/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                            W(i+1,NewCStateDec+1) = p15;
                            
                            %                             p17 = xi/Z;
                            p17 = lambda/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                            W(i+1,NewCStateDec+1) = p17;
                            
                        case '2'
                            
                            p24 = nu/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            W(i+1,NewCStateDec+1) = p24;
                            
                            %                             p27 = k/Z;
                            p27 = r/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                            W(i+1,NewCStateDec+1) = p27;
                            
                            p28 = lambda/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                            W(i+1,NewCStateDec+1) = p28;
                            
                        case '3'
                            
                            p34 = nu/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            W(i+1,NewCStateDec+1) = p34;
                            
                            %                             p36 = k/Z;
                            p36 = r/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                            W(i+1,NewCStateDec+1) = p36;
                            
                            p38 = lambda/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                            W(i+1,NewCStateDec+1) = p38;
                            
                        case '4'
                            
                            p40 = k/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                            W(i+1,NewCStateDec+1) = p40;
                            
                            p41 = k/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                            W(i+1,NewCStateDec+1) = p41;
                            
                            p42 = xi/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                            W(i+1,NewCStateDec+1) = p42;
                            
                            p43 = xi/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                            W(i+1,NewCStateDec+1) = p43;
                            
                    end
                end
                if(NumConnected)
                    for j=1:NumNonConnected
                        BCRj = TransitionNonConIndex(j);
                        BCRState = ClusterState(BCRj);
                        
                        switch BCRState
                            case '5' % transition to 5, 6, 7
                                p50 = q/Z;
                                NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                                W(i+1,NewCStateDec+1) = p50;
                                
                                p51 = q/Z;
                                NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                                W(i+1,NewCStateDec+1) = p51;
                            case '6'
                                p63 = q/Z;
                                NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                                W(i+1,NewCStateDec+1) = p63;
                            case '7'
                                p72 = q/Z;
                                NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                                W(i+1,NewCStateDec+1) = p72;
                        end
                    end
                end
                if( round(sum(W(i+1,:)),2)>1 )
                    sum(W(i+1,:))
                end
            end
            obj.W = W;
            obj.TCstates = TCstates;
            obj.ACstates = ACstates;
            
        end
        
        
        %%%%%
        %%%%%
        %%%%%
        
        function ConstructEQWMat(obj)
            %             NState = 2^(2*obj.NBCR);
            NAg = obj.NAg;
            kb = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            NBCR = obj.NBCR;
            %             F = 1e-12;
            %             F = obj.F;
            F=0;
            % Each BCR is a word (16 bit). The transition states are
            Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            StatesNum = length(Tstate)+length(Astate);
            obj.StatesNum = StatesNum;
            W = sparse(StatesNum^NBCR,StatesNum^NBCR);
            Wrates = sparse(StatesNum^NBCR,StatesNum^NBCR);
            EqState = [];
            ClusterEqStates = {};
            for i=0:(StatesNum^NBCR-1)
                ClusterState = dec2base(i,StatesNum,NBCR);
                
                TransitionIndex = regexp(ClusterState,'[0-4]');
                TransitionNonConIndex = regexp(ClusterState,'[5-7]');
                NumConnected = length(TransitionIndex);
                NumNonConnected = length(TransitionNonConIndex);
                if(NumConnected)
                    TCstates = [TCstates , i+1];
                else
                    ACstates = [ACstates , i+1];
                end
                %Connected states:
                f = F/NumConnected;
                
                lambda = exp(-obj.EAgMem)*exp(xb*f/(kb*T));
                r = exp(-obj.EAbAg)*exp(xb*f/(kb*T));
                k = exp(-obj.EAbAg)*exp(0.5*xb*f/(kb*T));
                xi =  exp(-obj.EAgMem)*exp(0.5*xb*f/(kb*T));
                
                nu = 0;
                
                NAgA = obj.NAgAttached(ClusterState);
                
                FreeAg = NAg-NAgA;
                if(FreeAg<0)
                    FreeAg=0; % cannot be in this state if NAg<NAgA
                    continue;
                end
                q = FreeAg*obj.q;
                %                 q = (NAg-1)*obj.q;
                %                 Z = obj.StatePartitionEq(ClusterState,f,FreeAg,NumConnected);
                Z=1;
                idxNoNoStaes = regexp(ClusterState,{'2' '3' '6' '7' '8'});
                if(find(~cellfun(@isempty,idxNoNoStaes)))
                    continue;
                end
                for j=1:NumConnected
                    BCRj = TransitionIndex(j);
                    BCRState = ClusterState(BCRj);
                    switch BCRState
                        case '0' % transition to 5, 6, 7
                            p04 = q/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            W(i+1,NewCStateDec+1) = p04;
                            
                            p05 = r/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                            W(i+1,NewCStateDec+1) = p05;
                            
                            %                             p06 = xi/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                            %                             W(i+1,NewCStateDec+1) = p06;
                        case '1'
                            p14 = q/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            W(i+1,NewCStateDec+1) = p14;
                            
                            p15 = r/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                            W(i+1,NewCStateDec+1) = p15;
                            %
                            %                             p17 = xi/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                            %                             W(i+1,NewCStateDec+1) = p17;
                            
                            %                         case '2'
                            %
                            %                             p24 = nu/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            %                             W(i+1,NewCStateDec+1) = p24;
                            %
                            %                             p27 = k/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                            %                             W(i+1,NewCStateDec+1) = p27;
                            %
                            %                             p28 = lambda/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                            %                             W(i+1,NewCStateDec+1) = p28;
                            %
                            %                         case '3'
                            %
                            %                             p34 = nu/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                            %                             W(i+1,NewCStateDec+1) = p34;
                            %
                            %                             p36 = k/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                            %                             W(i+1,NewCStateDec+1) = p36;
                            %
                            %                             p38 = lambda/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                            %                             W(i+1,NewCStateDec+1) = p38;
                            
                        case '4'
                            
                            p40 = k/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                            W(i+1,NewCStateDec+1) = p40;
                            
                            p41 = k/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                            W(i+1,NewCStateDec+1) = p41;
                            
                            %                             p42 = xi/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                            %                             W(i+1,NewCStateDec+1) = p42;
                            %
                            %                             p43 = xi/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                            %                             W(i+1,NewCStateDec+1) = p43;
                            
                    end
                end
                %                 if(NumConnected)
                for j=1:NumNonConnected
                    BCRj = TransitionNonConIndex(j);
                    BCRState = ClusterState(BCRj);
                    
                    switch BCRState
                        case '5' % transition to 5, 6, 7
                            p50 = q/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                            W(i+1,NewCStateDec+1) = p50;
                            
                            p51 = q/Z;
                            NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                            W(i+1,NewCStateDec+1) = p51;
                            %                         case '6'
                            %                             p63 = q/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                            %                             W(i+1,NewCStateDec+1) = p63;
                            %                         case '7'
                            %                             p72 = q/Z;
                            %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                            %                             W(i+1,NewCStateDec+1) = p72;
                    end
                end
                if(sum(W(i+1,:)))
                    EqState = [EqState ; i+1];
                    ClusterEqStates{length(ClusterEqStates)+1} = ClusterState;
                    %                     if( (i+1)==48)
                    %                         ClusterState
                    %                     end
                end
                
                Wrates(i+1,i+1) = -sum(W(i+1,:));
                idx = find(W(i+1,:));
                for ll=1:length(idx)
                    Wrates(idx(ll),i+1) = W(i+1,idx(ll));
                end
                W(i+1,:);
                %                 end
                sum(W(i+1,:));
                if( round(sum(W(i+1,:)),2)>1 )
                    sum(W(i+1,:));
                end
            end
            M = full(Wrates);
            K = M(EqState,EqState);
            [V D] = eig(K);
            %             K = [];
            %             for ll=1:length(M)
            %                 idx = find(M(ll,:));
            %
            %             end
            [v EqEigenidx] = min(diag(abs(D)));
            obj.Weq = K;
            
            StartingStateNames = {};
            StartingState = [];
            for i=1:length(ClusterEqStates)
                %                 NAgC = 0;
                %                 ClusterState = dec2base(TCstates(i)-1,obj.StatesNum,obj.NBCR);
                ClusterState = ClusterEqStates{i};
                BoundArms = regexp(ClusterState,'[0-1,4]');
                NonBoundAbs = regexp(ClusterState,'[5]');
                if( (length(BoundArms) > 0))
                    StartingState = [StartingState , EqState(i)];
                    c = length(StartingStateNames);
                    StartingStateNames{c+1} = ClusterEqStates{i};
                end
            end
            
            
            %             obj.StartingStateEq = EqState;
            %             obj.StartingStateNamesEq = ClusterEqStates;
            obj.StartingStateEq = StartingState;
            obj.StartingStateNamesEq = StartingStateNames;
            obj.Parr = V(:,EqEigenidx)/sum(V(:,EqEigenidx));
            %             obj.TCstates = TCstates;
            %             obj.ACstates = ACstates;
            
        end
        
        %%%%%
        %%%%%
        %%%%%
        
        function CaptureProb(obj)
            W = obj.W;
            TCstates = obj.TCstates;
            ACstates = obj.ACstates;
            tN = length(TCstates);
            rN = length(ACstates);
            
            Q = W(TCstates,TCstates);
            R = W(TCstates,ACstates);
            %             ZeroMat = zeros(rN,tN);
            %             I = eye(rN);
            
            It = speye(size(Q,1));
            NMat = inv(It-Q);
            B = NMat*R;
            
            
            ACstates_AgCap = obj.NAgCaptured(ACstates);
            AgCapPro = [];
            for i=1:size(B,1)
                for AgC=0:max(ACstates_AgCap)
                    idx = find(ACstates_AgCap==AgC);
                    AgCapPro(i,AgC+1) = sum(B(i,idx));
                end
            end
            obj.AgCapPro = AgCapPro;
        end
        
        
        function CaptureProb2Ag(obj)
            for n_bcr=1:obj.aff.NBCR
                W = obj.WnBCR{n_bcr};
                %             W = obj.W;
                TCstates = obj.TCstates;
                ACstates = obj.ACstates;
                tN = length(TCstates);
                rN = length(ACstates);
                
                Q = W(TCstates,TCstates);
                R = W(TCstates,ACstates);
                %             ZeroMat = zeros(rN,tN);
                %             I = eye(rN);
                
                It = speye(size(Q,1));
                NMat = inv(It-Q);
                B = NMat*R;
                
                
% %                 % The amount of capture Ag in each final state
%                 ACstates_AgCap = obj.NAgCaptured2Ag(ACstates);
%                 [C,ia,ic] = unique(ACstates_AgCap','rows');
%                 AgCapPro = zeros(size(B,1),length(C));
% %                 AgCapPro = [];
%                 [row,col] = find(B);
%                 rowB = unique(row);
%                 for i=rowB'
%                     for j=1:length(C)
%                         idx = find(ic == j);
%                         AgCapPro(i,j) = sum(B(i,idx))';
%                     end
%                 end
                % A shorter way to compute the amount of capture Ag in each final state
                if(obj.AgCapProbSumFlag==0)
                    obj.AgCapProbSumCalc(ACstates);
                end
%                 CalcAgCapProMat = B*obj.SumCapProbMat;
%                 AAA = AgCapPro-CalcAgCapProMat;
%                 if( length(find(AAA(:))) )
%                     AAA
%                 end
                AgCapPro = B*obj.SumCapProbMat;


%                 AgCapPro = [];
%                 for i=1:size(B,1) %03/12/18
%                     for j=1:length(C)
%                         idx = find(ic == j);
%                         AgCapPro(i,j) = sum(B(i,idx))';
%                     end
%                 end
                %             for i=1:size(B,1)
                %                 for AgC=0:max(ACstates_AgCap)
                %                     idx = find(ACstates_AgCap==AgC);
                %                     AgCapPro(i,AgC+1) = sum(B(i,idx));
                %                 end
                %             end
%                 The next line is only initialize once using the function AgCapProbSumCalc
%                 obj.AgCapState = C;  % the first column correspond to the amount of captured Ag of type 1 ; the second column correspond to the amount of Ag mol of type 2
                obj.AgCapProArr{n_bcr} = AgCapPro;
            end
        end
        
        function Z= StatePartition(obj,ClusterState,F,FreeAg,NumConnected)
            Z = 0;
            NAg = obj.NAg;
            kb = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            lambda = exp(-obj.EAgMem)*exp(xb*F/(kb*T));
            r = exp(-obj.EAbAg)*exp(xb*F/(kb*T));
            k = exp(-obj.EAbAg)*exp(0.5*xb*F/(kb*T));
            xi =  exp(-obj.EAgMem)*exp(0.5*xb*F/(kb*T));
            
            NAgA = obj.NAgAttached(ClusterState);
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
                    case '0'
                        %                         Z = Z + (r+q+xi);
                        Z = Z + (r+q+lambda);
                    case '1'
                        %                         Z = Z + (r+q+xi);
                        Z = Z + (r+q+lambda);
                    case '2'
                        %                         Z = Z + (k+nu+lambda);
                        Z = Z + (r+nu+lambda);
                    case '3'
                        %                         Z = Z + (k+nu+lambda);
                        Z = Z + (r+nu+lambda);
                    case '4'
                        Z = Z + (2*k+2*xi);
                    case '5'
                        if(NumConnected)
                            Z = Z + 2*q;
                        end
                    case '6'
                        if(NumConnected)
                            Z = Z + q;
                        end
                    case '7'
                        if(NumConnected)
                            Z = Z + q;
                        end
                end
            end
        end
        
        %%%%%
        %%%%%
        %%%%%
        
        function Z= StatePartitionEq(obj,ClusterState,F,FreeAg,NumConnected)
            Z = 0;
            NAg = obj.NAg;
            kb = 1.38064852e-23;
            T=273;
            xb = 1e-9;
            
            lambda = exp(-obj.EAgMem)*exp(xb*F/(kb*T));
            r = exp(-obj.EAbAg)*exp(xb*F/(kb*T));
            k = exp(-obj.EAbAg)*exp(0.5*xb*F/(kb*T));
            xi =  exp(-obj.EAgMem)*exp(0.5*xb*F/(kb*T));
            
            NAgA = obj.NAgAttached(ClusterState);
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
                    case '0'
                        %                         Z = Z + (r+q+xi);
                        Z = Z + (r+q);
                    case '1'
                        %                         Z = Z + (r+q+xi);
                        Z = Z + (r+q);
                        %                     case '2'
                        %                         Z = Z + (k+nu+lambda);
                        %                     case '3'
                        %                         Z = Z + (k+nu+lambda);
                    case '4'
                        %                         Z = Z + (2*k+2*xi);
                        Z = Z + (2*k);
                    case '5'
                        %                         if(NumConnected)
                        Z = Z + 2*q;
                        %                         end
                        %                     case '6'
                        % %                         if(NumConnected)
                        %                             Z = Z + q;
                        % %                         end
                        %                     case '7'
                        % %                         if(NumConnected)
                        %                             Z = Z + q;
                        % %                         end
                end
            end
        end
        
        
        %%%%%
        %%%%%
        %%%%%
        
        function NewCStateDec = ClusterStateNew(obj,ClusterState,idx,NewS);
            NewCState = ClusterState;
            NewCState(idx) = NewS;
            NewCStateDec = base2dec(NewCState,obj.StatesNum);
        end
        function NAgA = NAgAttached(obj,ClusterState)
            NAgA = 0;
            for n=1:obj.NBCR
                BCRState = ClusterState(n);
                switch BCRState
                    case '0'
                        NAgA = NAgA + 1;
                    case '1'
                        NAgA = NAgA + 1;
                    case '2'
                        NAgA = NAgA + 2;
                    case '3'
                        NAgA = NAgA + 2;
                    case '4'
                        NAgA = NAgA + 2;
                    case '6'
                        NAgA = NAgA + 1;
                    case '7'
                        NAgA = NAgA + 1;
                    case '8'
                        NAgA = NAgA + 2;
                end
            end
        end
        
        
        function [NAgaA NAgbA] = NAgAttached2Ag(obj,ClusterState)
            NAgaA = 0;
            NAgbA = 0;
            for n=1:obj.NBCR
                %                 BCRState = ClusterState(n);
                BCRState = ClusterState;
                switch BCRState
                    case '0'
                        NAgaA = NAgaA + 1;
                    case '1'
                        NAgbA = NAgbA + 1;
                    case '2'
                        NAgaA = NAgaA + 1;
                    case '3'
                        NAgaA = NAgaA + 2;
                    case '4'
                        NAgaA = NAgaA + 2;
                    case '5'
                        NAgaA = NAgaA + 2;
                    case '6'
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case '7'
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case '8'
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case '9'
                        NAgbA = NAgbA + 1;
                    case '10'
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case '11'
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case '12'
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case '13'
                        NAgbA = NAgbA + 2;
                    case '14'
                        NAgbA = NAgbA + 2;
                    case '15'
                        NAgbA = NAgbA + 2;
                end
            end
        end
        
%          for i=0:(StatesNum^NBCR-1)
%                     ClusterState = num2str(i);
        
        function [NAgaA NAgbA] = NAgAttached2AgInt(obj,ClusterState)
            NAgaA = 0;
            NAgbA = 0;
            for n=1:obj.NBCR
                %                 BCRState = ClusterState(n);
                BCRState = ClusterState;
                switch BCRState
                    case 0
                        NAgaA = NAgaA + 1;
                    case 1
                        NAgbA = NAgbA + 1;
                    case 2
                        NAgaA = NAgaA + 1;
                    case 3
                        NAgaA = NAgaA + 2;
                    case 4
                        NAgaA = NAgaA + 2;
                    case 5
                        NAgaA = NAgaA + 2;
                    case 6
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 7
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 8
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 9
                        NAgbA = NAgbA + 1;
                    case 10
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 11
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 12
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 13
                        NAgbA = NAgbA + 2;
                    case 14
                        NAgbA = NAgbA + 2;
                    case 15
                        NAgbA = NAgbA + 2;
                end
            end
        end
        
        function ACstates_AgCap = NAgCaptured(obj,ACstates)
            
            ACstates_AgCap = zeros(1,length(ACstates));
            for i=1:length(ACstates_AgCap)
                NAgC = 0;
                ClusterState = dec2base(ACstates(i)-1,obj.StatesNum,obj.NBCR);
                for n=1:size(ClusterState,2)
                    %                 ACstates,StatesNum,obj.NBCR)
                    BCRState = ClusterState(n);
                    switch BCRState
                        case '5'
                            NAgC = NAgC + 0;
                        case '6'
                            NAgC = NAgC + 1;
                        case '7'
                            NAgC = NAgC + 1;
                        case '8'
                            NAgC = NAgC + 2;
                    end
                end
                ACstates_AgCap(i) = NAgC;
            end
        end
        
        
        function ACstates_AgCap = NAgCaptured2Ag(obj,ACstates)
            %             Astate = [40,112,201,256,328,417,523,595,684];
            %             AstateName = {'001111', '011011', '021110', '100111', '110011', '120110', '201101', '211001', '221100'};
            
            ACstates_AgCap = zeros(2,length(ACstates));
            for i=1:length(ACstates_AgCap)
                %                 NAgC = 0;
                NAgaA = 0;
                NAgbA = 0;
  
%                 BCRState = num2str(ACstates(i));
                BCRState = ACstates(i);
                switch BCRState
                    case 17
                        0;
                    case 18
                        NAgaA = NAgaA + 1;
                    case 19
                        NAgbA = NAgbA + 1;
                    case 20
                        NAgaA = NAgaA + 1;
                    case 21
                        NAgaA = NAgaA + 2;
                    case 22
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 23
                        NAgbA = NAgbA + 1;
                    case 24
                        NAgaA = NAgaA + 1;
                        NAgbA = NAgbA + 1;
                    case 25
                        NAgbA = NAgbA + 2;
                        
                end
                %                 end
                %                 ACstates_AgCap(i) = NAgC;
                ACstates_AgCap(1,i) = NAgaA;
                ACstates_AgCap(2,i) = NAgbA;
            end
        end
        
        %%%%%%%
        
        function StartingState = StartingStatesFind2Ag(obj,TCstates)
            StartingState = [];
            StartingStateNames = {};
            %             ACstates_AgCap = zeros(1,length(ACstates));
            for i=1:length(TCstates)
                %                 NAgC = 0;
                %                 ClusterState = dec2base(TCstates(i)-1,obj.StatesNum,obj.NBCR);
                ClusterState = obj.TstateName{i}(1:2);
                BoundArms = regexp(ClusterState,'[1-2]');
                NonBoundAarms = regexp(ClusterState,'[0]');
                if( (length(BoundArms) == 1) & (length(NonBoundAarms) == 1))
                    StartingState = [StartingState , TCstates(i)];
                    c = length(StartingStateNames);
                    StartingStateNames{c+1} = obj.TstateName{i};
                end
                %                 BoundArms = regexp(ClusterState,'[0-1]');
                %                 NonBoundAbs = regexp(ClusterState,'[5]');
                %                 if( (length(BoundArms) == 1) & (length(NonBoundAbs) == obj.NBCR-1))
                %                     StartingState = [StartingState , TCstates(i)];
                %                     c = length(StartingStateNames);
                %                     StartingStateNames{c+1} = ClusterState;
                %                 end
            end
            obj.StartingState = StartingState;
            obj.StartingStateNames = StartingStateNames;
        end
        
        %%%%
        
        function StartingState = StartingStatesFind(obj,TCstates)
            StartingState = [];
            StartingStateNames = {};
            %             ACstates_AgCap = zeros(1,length(ACstates));
            for i=1:length(TCstates)
                %                 NAgC = 0;
                ClusterState = dec2base(TCstates(i)-1,obj.StatesNum,obj.NBCR);
                BoundArms = regexp(ClusterState,'[0-1]');
                NonBoundAbs = regexp(ClusterState,'[5]');
                if( (length(BoundArms) == 1) & (length(NonBoundAbs) == obj.NBCR-1))
                    StartingState = [StartingState , TCstates(i)];
                    c = length(StartingStateNames);
                    StartingStateNames{c+1} = ClusterState;
                end
            end
            obj.StartingState = StartingState;
            obj.StartingStateNames = StartingStateNames;
        end
        
        
        %%%%%
        %%%%%
        
        function AbsorbingStatesFind(obj)
            
            Tstate = [7,11,13,14,15]; % 0111, 1011, 1101, 1110, 1111 % connected states : force falls on the bond
            Astate = [3,6,9,12]; % 0011 Rap, 0110 Cap 1 arm, 1001 Cap 1 arm, 1100 Cap 2 arms
            TCstates = []; % Transition cluster states
            ACstates = []; % Absorbing cluster states
            StatesNum = length(Tstate)+length(Astate);
            obj.StatesNum = StatesNum;
            %             W = sparse(StatesNum^NBCR,StatesNum^NBCR);
            NBCR = obj.NBCR;
            
            for i=0:(StatesNum^NBCR-1)
                ClusterState = dec2base(i,StatesNum,NBCR);
                
                TransitionIndex = regexp(ClusterState,'[0-4]');
                TransitionNonConIndex = regexp(ClusterState,'[5-7]');
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
        
        function CalcProbDensityNoneq2Ag(obj);
            
            global ICArr;
                        
            IC = ICArr{obj.ICnum(1)}{obj.ICnum(2)};
            
            for n_bcr=1:obj.aff.NBCR
            Parr = [];
%             NAga = IC(1);
%             NAgb = IC(2);
            NAga = IC{n_bcr}(1);
            NAgb = IC{n_bcr}(2);
            n=NAga+NAgb;
            %             w = obj.w;
            tau = obj.tau;
            
            %             n=obj.NAg;
            
            % We start from 00 and have a chance of getting to 10/01/11
            % before time tau;
            %             qa = obj.qa;
            %             qb = obj.qb;
            qa=1;
            qb=1;
            switch n
                case 0
                    Parr(1) = 1;
                    Parr(2) = 0;
                    Parr(3) = 0;
                    Parr(4) = 0;
                    Parr(5) = 0;
                otherwise
                    Parr(1) = exp(- 2*tau*(qa*obj.NBCR*NAga + qb*obj.NBCR*NAgb) );
                    
                    %                     Parr(1) = exp(-2*q*obj.NBCR*obj.NAg*tau);
                    POneArm = 1-Parr(1);
                    %                     p=POneArm/(2*obj.NBCR);
                    Z = 2*(qa*NAga+qb*NAgb);
                    for j=1:length(obj.StartingStateNames)
                        BCRs = obj.StartingStateNames{j}(1:2);
                        BoundArm = regexp(BCRs,'[1-2]');
                        if(BoundArm)
                            if(BCRs(BoundArm)=='1')
                                %                         Parr(j+1) = p;
                                Parr(j+1) = POneArm*qa*NAga/Z;
                            else
                                Parr(j+1) = POneArm*qb*NAgb/Z;
                            end
                        end
                    end
                    %                     Parr(2) = 0.5*(1-exp(-2*q*n*tau));
                    %                     Parr(3) = 0.5*(1-exp(-2*q*n*tau));
                    %                     Parr(4) = 0;
            end
            obj.ParrNBCR{n_bcr} = Parr;
            end
        end
        
        %%%%
        %%%%
        %%%%
        %%%%
        %%%%
        %%%%
        
        function CalcProbDensityNoneq2AgExplicitDensity(obj);
            
            global ICArr;
                        
            IC = ICArr{obj.ICnum(1)}{obj.ICnum(2)};
            
            for n_bcr=1:obj.aff.NBCR
            Parr = [];
%             NAga = IC(1);
%             NAgb = IC(2);
            NAga = IC{n_bcr}(1);
            NAgb = IC{n_bcr}(2);
            n=NAga+NAgb;
            %             w = obj.w;
            tau = obj.tau;
            
            % We start from 00 and have a chance of getting to 10/01/11
            % before time tau;
            %             qa = obj.qa;
            %             qb = obj.qb;            
            
            qa=1;
            qb=1;
%             Epitope = 1;
%             qaNAga = DensityDependentRate(NAga,obj.aff.Epitope,1); % The on rate for the first arm to antigen a
%             qbNAgb = DensityDependentRate(NAgb,obj.aff.Epitope,1); % The on rate for the first arm to antigen b
            
            
            qaNAga = obj.DensityDependentRate1Arm(NAga,obj.aff.Epitope); % binding rate to Ag "a" for the first arm
            qbNAgb = obj.DensityDependentRate1Arm(NAgb,obj.aff.Epitope); % binding rate to Ag "b" for the first arm

            
            switch n
                case 0
                    Parr(1) = 1;
                    Parr(2) = 0;
                    Parr(3) = 0;
                    Parr(4) = 0;
                    Parr(5) = 0;
                otherwise
%                     Parr(1) = exp(- 2*tau*(qa*obj.NBCR*NAga + qb*obj.NBCR*NAgb) );
                    Parr(1) = exp(- tau*(qaNAga + qbNAgb) );

                    POneArm = 1-Parr(1);

%                     Z = 2*(qa*NAga+qb*NAgb);
                    Z = qaNAga + qbNAgb;
                    for j=1:length(obj.StartingStateNames)
                        BCRs = obj.StartingStateNames{j}(1:2);
                        BoundArm = regexp(BCRs,'[1-2]');
                        if(BoundArm)
                            if(BCRs(BoundArm)=='1')
                                %                         Parr(j+1) = p;
%                                 Parr(j+1) = POneArm*qa*NAga/Z;
                                Parr(j+1) = POneArm*0.5*qaNAga/Z;
                            else
%                                 Parr(j+1) = POneArm*qb*NAgb/Z;
                                Parr(j+1) = POneArm*0.5*qbNAgb/Z;
                            end
                        end
                    end
            end
            obj.ParrNBCR{n_bcr} = Parr;
            end
        end
        
        %%%%
        %%%%
        %%%%
        %%%%
        %%%%
        %%%%
        
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
        function [NextStateName NextStateNum NextStateProb] = TransitionSimulate(obj,ClusterState);
            
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
            
            TransitionIndex = regexp(ClusterState,'[0-4]');
            TransitionNonConIndex = regexp(ClusterState,'[5-7]');
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
            
            NAgA = obj.NAgAttached(ClusterState);
            
            FreeAg = NAg-NAgA;
            if(FreeAg<0)
                FreeAg=0; % cannot be in this state if NAg<NAgA
                %                 continue;
            end
            q = FreeAg*obj.q;
            Z = obj.StatePartition(ClusterState,f,FreeAg,NumConnected);
            
            NextStateNum = [];
            NextStateProb = [];
            NextStateName = {};
            
            for j=1:NumConnected
                BCRj = TransitionIndex(j);
                BCRState = ClusterState(BCRj);
                switch BCRState
                    case '0' % transition to 5, 6, 7
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'5',BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'6',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        %                             p04 = q/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                        %                             W(i+1,NewCStateDec+1) = p04;
                        %
                        %                             p05 = r/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                        %                             W(i+1,NewCStateDec+1) = p05;
                        %
                        %                             p06 = xi/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                        %                             W(i+1,NewCStateDec+1) = p06;
                    case '1'
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'5',BCRj,r/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'7',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        %                             p14 = q/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                        %                             W(i+1,NewCStateDec+1) = p14;
                        %
                        %                             p15 = r/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'5');
                        %                             W(i+1,NewCStateDec+1) = p15;
                        %
                        %                             p17 = xi/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                        %                             W(i+1,NewCStateDec+1) = p17;
                        
                    case '2'
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,nu/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'7',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'8',BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);
                        %                             p24 = nu/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                        %                             W(i+1,NewCStateDec+1) = p24;
                        %
                        %                             p27 = k/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'7');
                        %                             W(i+1,NewCStateDec+1) = p27;
                        %
                        %                             p28 = lambda/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                        %                             W(i+1,NewCStateDec+1) = p28;
                        
                    case '3'
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'4',BCRj,nu/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'6',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'8',BCRj,lambda/Z,NextStateNum,NextStateProb,NextStateName);
                        %                             p34 = nu/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'4');
                        %                             W(i+1,NewCStateDec+1) = p34;
                        %
                        %                             p36 = k/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'6');
                        %                             W(i+1,NewCStateDec+1) = p36;
                        %
                        %                             p38 = lambda/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'8');
                        %                             W(i+1,NewCStateDec+1) = p38;
                        
                    case '4'
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'0',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'1',BCRj,k/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'2',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'3',BCRj,xi/Z,NextStateNum,NextStateProb,NextStateName);
                        
                        %                             p40 = k/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                        %                             W(i+1,NewCStateDec+1) = p40;
                        %
                        %                             p41 = k/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                        %                             W(i+1,NewCStateDec+1) = p41;
                        %
                        %                             p42 = xi/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                        %                             W(i+1,NewCStateDec+1) = p42;
                        %
                        %                             p43 = xi/Z;
                        %                             NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                        %                             W(i+1,NewCStateDec+1) = p43;
                        
                end
            end
            if(NumConnected)
                for j=1:NumNonConnected
                    BCRj = TransitionNonConIndex(j);
                    BCRState = ClusterState(BCRj);
                    
                    switch BCRState
                        case '5' % transition to 5, 6, 7
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'0',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'1',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            
                            %                                 p50 = q/Z;
                            %                                 NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'0');
                            %                                 W(i+1,NewCStateDec+1) = p50;
                            %
                            %                                 p51 = q/Z;
                            %                                 NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'1');
                            %                                 W(i+1,NewCStateDec+1) = p51;
                        case '6'
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'3',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            
                            %                                 p63 = q/Z;
                            %                                 NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'3');
                            %                                 W(i+1,NewCStateDec+1) = p63;
                        case '7'
                            [NextStateName NextStateNum NextStateProb] = obj.UpdateNextState(ClusterState,'2',BCRj,q/Z,NextStateNum,NextStateProb,NextStateName);
                            %                                 p72 = q/Z;
                            %                                 NewCStateDec = obj.ClusterStateNew(ClusterState,BCRj,'2');
                            %                                 W(i+1,NewCStateDec+1) = p72;
                    end
                end
            end
            if( round(sum(NextStateProb),2)>1 )
                sum(NextStateProb)
            end
            
            %             obj.W = W;
            %             obj.TCstates = TCstates;
            %             obj.ACstates = ACstates;
        end
        
        function [EndStateName EndStateNum ACstates_AgCap]  = SimulateClusterAbs(obj,InitialState)
            
            NextStateNum = base2dec(InitialState,obj.StatesNum);
            NextStateName = InitialState;
            count = 0;
            while ( length(find(obj.ACstates == NextStateNum+1)) == 0)
                [NextStateNameArr NextStateNumArr NextStateProbArr] = obj.TransitionSimulate(NextStateName);
                
                P = cumsum(NextStateProbArr);
                ii = min(find(P-rand>0));
                NextStateName = NextStateNameArr{ii};
                NextStateNum = NextStateNumArr(ii);
                count = count+1;
            end
            EndStateName = NextStateName;
            EndStateNum = NextStateNum;
            ACstates_AgCap = NAgCaptured(obj,EndStateNum+1);
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
        function AgCapProbSumCalc(obj,ACstates)
            ACstates_AgCap = obj.NAgCaptured2Ag(ACstates);
            AbsorbingStateNum = length(ACstates);
            [C,ia,ic] = unique(ACstates_AgCap','rows');
            AgCapResultsNum = length(C);
            
            SumCapProbMat = zeros(AbsorbingStateNum,AgCapResultsNum);
            for j=1:length(C)
                idx = find(ic == j);
                SumCapProbMat(idx,j) = 1;
            end
            obj.SumCapProbMat = SumCapProbMat;
            obj.AgCapState = C;
            obj.AgCapProbSumFlag = 1;
        end
        
        function kon = DensityDependentRate1Arm(obj,N,Epitope)
            konMean = obj.aff.KonArr.kOn1Arm(Epitope,N+1);
            p = obj.aff.KonArr.P1Arm(Epitope,N+1);
            r = rand;
            if(r<p)
                kon = konMean;
            else
                kon = 0;
            end
        end
        
        function kon = DensityDependentRate2Arm(obj,N,Epitope)
            konMean = obj.aff.KonArr.kOn2Arm(Epitope,N+1);
            p = obj.aff.KonArr.P2Arm(Epitope,N+1);
            r = rand;
            if(r<p)
                kon = konMean;
            else
                kon = 0;
            end
        end
        
    end
    
end
% end