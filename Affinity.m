classdef Affinity < handle
    properties
        NBCR; % number of antigens
        k;
        q;
        q0; %kon density independent
        q1;
        q2;
        w;
        w1;
        w2;
        BondsNum;
        Parr;
        EAbAg;
        EAgMem;
        SecBind = 2; % effect of cooperativiate
        F;
        a;
        A;
        NAg;
        BindFlag;
        tau;
        AbAgClusterProp;
        CapSimFlag;
        EqFlag;
        MAg;
        SynapseRandomFlag;
        %         paramFlag;
        aff;
    end
    methods
        function obj = Affinity(v,aff);
            if(nargin > 0) %paramFlag, specify which is the input parameter: 0 w (affinity) / 1 q (kon) / 2 k (koff) / 3 E0 interaction energy
                
                %                 switch aff.paramFlag
                %                     case 0
                %                         obj.w = v;
                %                         obj.q = v;
                %                         obj.EAbAg = aff.EAbAg;
                %                         obj.k = obj.q/obj.w;
                %                     case 1
                %                         obj.w = v;
                %                         obj.q = v;
                %                         obj.k = obj.q/obj.w;
                %                     case 2
                %                         obj.q = 1;
                %                         obj.k = v;
                %                         obj.w = obj.q/obj.k;
                %                     case 3
                %                         %                         obj.q = 5;
                %                         obj.q = aff.q;
                %                         obj.EAbAg = v;
                %                         obj.k = exp(-v);
                %                         obj.w = obj.q/obj.k;
                %                     case 4
                %                         obj.NAg = aff.NAg;
                %                         obj.q = obj.NAg*aff.q;
                %                         obj.EAbAg = v;
                %                         obj.k = exp(-v);
                %                         obj.w = obj.q/obj.k;
                %                     case 5
                %                         obj.NAg = aff.NAg;
                %                         obj.q1 = obj.NAg*aff.q1;
                %                         obj.q2 = obj.NAg*v;
                %                         %                         obj.q2 = obj.NAg*aff.q2;
                %                         obj.EAbAg = aff.EAbAg;
                %                         obj.k = exp(-obj.EAbAg);
                %                         obj.w1 = obj.q1/obj.k;
                %                         obj.w2 = obj.q2/obj.k;
                %                     case 6
                %                         obj.NAg = aff.NAg;
                %                         obj.q0 = v;
                %                         %                         obj.q = obj.NAg*v;
                %                         obj.q = v;
                %                         obj.EAbAg = aff.EAbAg;
                %                         obj.k = exp(-obj.EAbAg);
                %                         obj.w = obj.q/obj.k;
                %                      case 7
                %                         obj.NAg = aff.NAg;
                %                         obj.q0 = aff.q1;
                %                         obj.q = aff.q1;
                %                         obj.NAg = aff.NAg;
                %                         obj.EAbAg = v;
                %                         obj.k = exp(-obj.EAbAg);
                %                         obj.w = obj.q/obj.k;
                %                     case 8
                %                         obj.NAg = aff.NAg;
                %                         obj.q1 = aff.q1;
                % %                         obj.q2 = aff.q2;
                %                         obj.q2 = v;
                %                         obj.q = obj.q2;
                %                         obj.NAg = aff.NAg;
                %                         obj.EAbAg = aff.EAbAg;
                %                         obj.k = exp(-obj.EAbAg);
                % %                         obj.w = obj.q/obj.k;
                %                         obj.w1 = obj.q1/obj.k;
                %                         obj.w2 = obj.q2/obj.k;
                %                     case 9
                %                         obj.NAg = aff.NAg;
                %                         obj.q1 = aff.q1;
                %                         obj.q2 = aff.q2;
                %                         obj.NAg = aff.NAg;
                %                         obj.EAbAg = aff.EAbAg;
                %                         obj.k = exp(-obj.EAbAg);
                %                         obj.w1 = obj.q1/obj.k;
                %                         obj.w2 = obj.q2/obj.k;
                %                 end
                
                obj.CapSimFlag = 0;
                obj.EqFlag = 0;
                obj.SynapseRandomFlag = 0; % When equal to zero, NAg is randomize once per B cell. When equal to 1, NAg is randomize for each cluster.
                
                obj.aff = aff;
                obj.NAg = aff.NAg;
                obj.q1 = aff.q1;
                obj.q2 = aff.q2;
                obj.NAg = aff.NAg;
                obj.EAbAg = aff.EAbAg;
                obj.k = exp(-obj.EAbAg);
                obj.w1 = obj.q1/obj.k;
                %                 obj.w2 = obj.q2/obj.k;
                obj.F = aff.F;
                
                obj.BindFlag = aff.BindFlag;
                obj.NBCR=aff.NBCR;
                obj.SecBind = aff.SecBind;
                obj.EAgMem = aff.EAgMem;
                obj.BondsNum=aff.BondsNum;
                obj.tau=aff.tau;
                %                 switch aff.paramFlag
                switch aff.BindFlag
                    case 1 % only 1 k_on ; two Ags
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                    case 2 % two k_on's ; two Ags
                        obj.q = aff.q2;
                    case 3 % one k_on ; one Ag
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                    case 4 % two k_on's ; one Ag
                        obj.q = aff.q2;
                    case 5 % Initial state 00;  only 1 k_on ; two Ags
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                    case 6 % Initial state 00;  only 1 k_on ; one Ag
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                    case 7 % Initial state 00; Multiple Ags ; BCR cluster
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                        if(obj.CapSimFlag)
                            obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize);
                        else
                            obj.AbAgClusterProp = AbAgCluster(obj.NBCR,obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize,obj.EqFlag);
                        end
                    case 8 % Initial state 00; Multiple Ags ; BCR cluster ; Ag distribution
                        obj.A = aff.A;
                        obj.AgDistribution(aff.NAg);
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                        if(obj.CapSimFlag)
                            obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize);
                        else
                            obj.AbAgClusterProp = AbAgCluster(aff,obj.NBCR,obj.NAg,obj.q,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize,obj.EqFlag);
                        end
                    case 9 % Initial state 00; Multiple Ags ; BCR cluster ; Ag distribution
                        obj.A = aff.A;
                        obj.AgDistribution(aff.NAg);
                        obj.q = aff.q1;
                        obj.w = obj.q/obj.k;
                        if(obj.CapSimFlag)
                            obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize);
                        else
                            obj.AbAgClusterProp = AbAgCluster(aff,obj.NBCR,obj.NAg,obj.q,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize,obj.EqFlag);
                        end
                        
                    case 10 % Initial state 00; Multiple Ags ; BCR cluster ; Ag distribution
                        % The approching time to the immune complex is a
                        % different from the one rate from the first or
                        % second BCR arm.
                        obj.A = aff.A;
                        %                         obj.AgDistribution(aff.NAg);
                        obj.AgDistributionNanoParticle(aff.NAg);
                        %                         qMFPTImmuneComp; % This is the arrival rate to
                        %                         the immune complex
                        obj.q = aff.q1;
                        obj.q1 = aff.q1;
                        obj.q2 = aff.q2;
                        obj.w = obj.q/obj.k;
                        if(obj.CapSimFlag)
                            obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize);
                        else
                            obj.AbAgClusterProp = AbAgCluster(aff,obj.NBCR,obj.NAg,obj.q,obj.EAbAg,obj.tau,obj.F,aff.ClusterSize,obj.EqFlag);
                        end
                end
                
                %                 switch aff.paramFlag
                switch aff.BindFlag
                    case 1
                        obj.CalcProbDensity();
                    case 2
                        obj.CalcProbDensity_q1q2();
                    case 3
                        obj.CalcProbDensity_SameAg();
                    case 4
                        obj.CalcProbDensity_q1q2SameAg();
                    case 5
                        obj.CalcProbDensityNoneq();
                    case 6
                        obj.CalcProbDensityNoneq();
                    case 7
                        obj.Parr = obj.AbAgClusterProp.Parr;
                    case 8
                        obj.Parr = obj.AbAgClusterProp.Parr;
                    case 9
                        obj.Parr = obj.AbAgClusterProp.Parr;
                    case 10
                        obj.Parr = obj.AbAgClusterProp.Parr;
                end
                
                %
                %                 switch obj.BindFlag
                %                     case 0
                %                         obj.CalcProb();
                %                     case 1
                %                         obj.CalcProbDensity();
                %                     case 2
                %                         obj.CalcProbDensity_q1q2();
                %                     case 3
                %                         obj.CalcProbDensity();
                %                 end
            else
            end
        end
        
        function CalcProbDensity(obj);
            
            w = obj.w;
            
            n=obj.NAg;
            
            switch n
                case 0
                    obj.Parr(1) = 1;
                    obj.Parr(2) = 0;
                    obj.Parr(3) = 0;
                    obj.Parr(4) = 0;
                case 1
                    N = 2+1/w;
                    obj.Parr(1) = 1/w;
                    obj.Parr(2) = 1;
                    obj.Parr(3) = 1;
                    obj.Parr(4) = 0;
                    obj.Parr = obj.Parr/N;
                otherwise
                    N = (1 + n*w*(2 + (n-1)*w))/((n-1)*n*w^2);
                    obj.Parr(1) = (n*(n-1)*w^2)^(-1);
                    obj.Parr(2) = ((n-1)*w)^(-1);
                    obj.Parr(3) = ((n-1)*w)^(-1);
                    obj.Parr(4) = 1;
                    obj.Parr = obj.Parr/N;
            end
        end
        
        %%%
        %%%
        
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
                    %                 case 1
                    %                     N = 2+1/w;
                    %                     obj.Parr(1) = 1/w;
                    %                     obj.Parr(2) = 1;
                    %                     obj.Parr(3) = 1;
                    %                     obj.Parr(4) = 0;
                    %                     obj.Parr = obj.Parr/N;
                otherwise
                    obj.Parr(1) = exp(-2*q*n*tau);
                    obj.Parr(2) = 0.5*(1-exp(-2*q*n*tau));
                    obj.Parr(3) = 0.5*(1-exp(-2*q*n*tau));
                    obj.Parr(4) = 0;
                    %                     N = (1 + n*w*(2 + (n-1)*w))/((n-1)*n*w^2);
                    %                     obj.Parr(1) = (n*(n-1)*w^2)^(-1);
                    %                     obj.Parr(2) = ((n-1)*w)^(-1);
                    %                     obj.Parr(3) = ((n-1)*w)^(-1);
                    %                     obj.Parr(4) = 1;
                    %                     obj.Parr = obj.Parr/N;
            end
        end
        
        %%%
        %%%
        
        function CalcProbDensity_SameAg(obj);
            
            %             w1 = obj.w1;
            %             w2 = obj.w2;
            w = obj.w;
            n=obj.NAg;
            
            switch n
                case 0
                    obj.Parr(1) = 1;
                    obj.Parr(2) = 0;
                    obj.Parr(3) = 0;
                    obj.Parr(4) = 0;
                otherwise
                    N = 1+ 2*(w)^-1 + (n*w^2)^(-1);
                    obj.Parr(1) = ( n*w^2 )^(-1);
                    obj.Parr(2) = (w)^(-1);
                    obj.Parr(3) = (w)^(-1);
                    obj.Parr(4) = 1;
                    obj.Parr = obj.Parr/N;
            end
        end
        
        %%%
        
        
        function CalcProbDensity_q1q2(obj);
            
            w1 = obj.w1;
            w2 = obj.w2;
            n=obj.NAg;
            
            switch n
                case 0
                    obj.Parr(1) = 1;
                    obj.Parr(2) = 0;
                    obj.Parr(3) = 0;
                    obj.Parr(4) = 0;
                case 1
                    N = 2+1/w1;
                    obj.Parr(1) = 1/w1;
                    obj.Parr(2) = 1;
                    obj.Parr(3) = 1;
                    obj.Parr(4) = 0;
                    obj.Parr = obj.Parr/N;
                otherwise
                    N = 1+ 2*((n-1)*w2)^-1 + (n*(n-1)*w1*w2)^(-1);
                    obj.Parr(1) = ( n*(n-1)*w1*w2 )^(-1);
                    obj.Parr(2) = ((n-1)*w2)^(-1);
                    obj.Parr(3) = ((n-1)*w2)^(-1);
                    obj.Parr(4) = 1;
                    obj.Parr = obj.Parr/N;
            end
        end
        
        %%%
        %%%
        
        function CalcProbDensity_q1q2SameAg(obj);
            
            w1 = obj.w1;
            w2 = obj.w2;
            n=obj.NAg;
            
            switch n
                case 0
                    obj.Parr(1) = 1;
                    obj.Parr(2) = 0;
                    obj.Parr(3) = 0;
                    obj.Parr(4) = 0;
                otherwise
                    N = 1+ 2*(w2)^-1 + (n*w1*w2)^(-1);
                    obj.Parr(1) = ( n*w1*w2 )^(-1);
                    obj.Parr(2) = (w2)^(-1);
                    obj.Parr(3) = (w2)^(-1);
                    obj.Parr(4) = 1;
                    obj.Parr = obj.Parr/N;
            end
        end
        
        %%%
        
        
        function CalcProb(obj);
            w = obj.w;
            switch obj.BondsNum
                case 1
                    obj.Parr(1) = (1+w)^(-1);
                    obj.Parr(2) = (1+1/w)^(-1);
                case 2
                    obj.Parr(1) = (w^2+2*w+1)^(-1);
                    obj.Parr(2) = (w+2+1/w)^(-1);
                    obj.Parr(3) = (w+2+1/w)^(-1);
                    obj.Parr(4) = (1+2/w+w^(-2))^(-1);
                    
            end
        end
        
        %%%
        function TotAg = CAgCaptureFunc(obj);
            
            switch obj.BindFlag
                case 1 % only 1 k_on
                    %                         TotAg = obj.CAgCaptureDensityCooperativity();
                    TotAg = obj.CAgCaptureDensity2Arms2Ags();
                case 2 % two k_on's
                    TotAg = obj.CAgCaptureDensity2Arms2Ags();
                case 3 % only 1 k_on
                    TotAg = obj.CAgCaptureDensityCooperativity();
                case 4 % two k_on's
                    TotAg = obj.CAgCaptureDensityCooperativity();
                case 5 % only 1 k_on
                    TotAg = obj.CAgCaptureDensity2Arms2Ags();
                case 6 % only 1 k_on Initial state 00 ; Same Ag
                    TotAg = obj.CAgCaptureDensityCooperativity();
                case 7 % Initial state 00; Multiple Ags ; BCR cluster
                    TotAg = obj.CAgCaptureDensityCluster();
                case 8 % Initial state 00; Multiple Ags ; BCR cluster; Ag distribution
                    obj.AgDistribution(obj.MAg);
                    %                     obj.q = aff.q1;
                    %                     obj.w = obj.q/obj.k;
                    if(obj.CapSimFlag)
                        obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,obj.AbAgClusterProp.ClusterSize);
                    else
                        obj.AbAgClusterProp = AbAgCluster(obj.aff,obj.NBCR,obj.NAg,obj.q,obj.EAbAg,obj.tau,obj.F,obj.AbAgClusterProp.ClusterSize,obj.EqFlag);
                    end
                    obj.Parr = obj.AbAgClusterProp.Parr;
                    if(obj.SynapseRandomFlag==0)
                        TotAg = obj.CAgCaptureDensityCluster();
                    else
                        TotAg = obj.CAgCaptureDensityClusterDistribution();
                    end
                case {9 , 10} % Initial state 00; Multiple Ags ; BCR cluster; Ag distribution; two Ag types.
                    %                     obj.AgDistribution(obj.MAg);
                    %                     obj.q = aff.q1;
                    %                     obj.w = obj.q/obj.k;
                    % Do we need to calculate this again????
                    if(obj.CapSimFlag)
                        obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,obj.AbAgClusterProp.ClusterSize);
                    else
                        obj.AbAgClusterProp = AbAgCluster(obj.aff,obj.NBCR,obj.NAg,obj.q,obj.EAbAg,obj.tau,obj.F,obj.AbAgClusterProp.ClusterSize,obj.EqFlag);
                    end
                    %                     obj.Parr = obj.AbAgClusterProp.Parr;
                    
                    if(obj.SynapseRandomFlag==0)
                        [TotAg1 TotAg2]  = obj.CAgCaptureDensityCluster2Ag();
                        TotAg = TotAg1+TotAg2;
                    else
                        TotAg = obj.CAgCaptureDensityClusterDistribution();
                    end
                    
            end
            
            %             switch obj.BindFlag
            %                 case 0
            %                     TotAg = obj.CAgCaptureDoubleBind();
            %                 case 1
            %                     TotAg = obj.CAgCaptureDensity();
            %                 case 2
            %                     TotAg = obj.CAgCaptureDensityCooperativity();
            %                 case 3
            %                     TotAg = obj.CAgCaptureDensityCooperativity();
            %             end
        end
        %%%
        function TotAg = CAgCaptureDoubleBind(obj);
            NBCR = obj.NBCR;
            TotAg = 0;
            P = cumsum(obj.Parr);
            for i=1:NBCR
                ii = min(find(P-rand>0));
                idx = strfind(dec2bin(ii-1,obj.BondsNum),'1');
                BitOn = length(idx);
                switch BitOn
                    case 0
                        E=0;
                    case 1
                        E = obj.EAbAg;
                    case 2
                        E = obj.SecBind*obj.EAbAg;
                    otherwise
                        E=0;
                end
                Z = exp(-obj.EAgMem) + exp(-E);
                if(E==0)
                    pRapAbAg = 0;
                else
                    pRapAbAg = exp(-obj.EAgMem)/Z;
                end
                if(rand<pRapAbAg)
                    TotAg = TotAg+1;
                end
            end
        end
        
        %%%
        
        function TotAg = CAgCaptureDensity(obj);
            NBCR = obj.NBCR;
            TotAg = 0;
            P = cumsum(obj.Parr);
            for i=1:NBCR
                ii = min(find(P-rand>0));
                idx = strfind(dec2bin(ii-1,obj.BondsNum),'1');
                BitOn = length(idx);
                switch BitOn
                    case 0
                        E=0;
                    case 1
                        E = obj.EAbAg;
                    case 2
                        %                         E = obj.SecBind*obj.EAbAg;
                        E = obj.EAbAg;
                    otherwise
                        E=0;
                end
                Z = exp(-obj.EAgMem) + exp(-E);
                if(E==0)
                    pRapAbAg = 0;
                else
                    pRapAbAg = exp(-obj.EAgMem)/Z;
                end
                for i=1:BitOn
                    if(rand<pRapAbAg)
                        TotAg = TotAg+1;
                    end
                end
            end
        end
        %%%
        
        
        
        %%%%%%%%%%%%
        %%%%%%%%%%%%
        
        function TotAg = CAgCaptureDensityCooperativity(obj);
            
            kb = 1.38064852e-23;
            %             F = 20e-12;
            %             F = 1e-12;
            F = obj.F;
            T=273;
            xb = 1e-9;
            xb*F/(kb*T);
            
            lambda = exp(-obj.EAgMem)*exp(xb*F/(kb*T));
            r = exp(-obj.EAbAg)*exp(xb*F/(kb*T));
            k = exp(-obj.EAbAg)*exp(0.5*xb*F/(kb*T));
            
            %             q = (obj.NAg-1)*obj.q;
            q = obj.q;
            %             switch obj.BindFlag
            %                 case 1
            %                     q = (obj.NAg-1)*obj.q;
            %                 case 2
            %                     q = (obj.NAg-1)*obj.q;
            %                 case 3
            %                     q = obj.q;
            %             end
            
            Cap7 = lambda.*(2.*k+q+r+lambda).*(2.*k.*(r+lambda)+lambda.*(q+r+lambda)).^(-1);
            Rap7 = 2.*k.*r.*(2.*k.*(r+lambda)+lambda.*(q+r+lambda)).^(-1);
            Cap5 = lambda.*(2.*k+q+lambda).*(2.*k.*(r+lambda)+lambda.*(q+r+lambda)).^(-1);
            Rap5 = r.*(2.*k+lambda).*(2.*k.*(r+lambda)+lambda.*(q+r+lambda)).^(-1);
            Cap3 = lambda.*(2.*k+q+lambda).*(2.*k.*(r+lambda)+lambda.*(q+r+lambda)).^(-1);
            Rap3 = r.*(2.*k+lambda).*(2.*k.*(r+lambda)+lambda.*(q+r+lambda)).^(-1);
            
            NBCR = obj.NBCR;
            TotAg = 0;
            P = cumsum(obj.Parr);
            for i=1:NBCR
                ii = min(find(P-rand>0));
                idx = strfind(dec2bin(ii-1,obj.BondsNum),'1');
                BitOn = length(idx);
                switch BitOn
                    case 0
                        %                         E=0;
                        Pcapture = 0;
                    case 1
                        %                         E = obj.EAbAg;
                        Pcapture = Cap5;
                    case 2
                        %                         E = obj.EAbAg;
                        Pcapture = Cap7;
                    otherwise
                        %                         E=0;
                        Pcapture = 1;
                end
                
                %                 Z = exp(-obj.EAgMem) + exp(-E);
                %                 if(E==0)
                %                     pRapAbAg = 0;
                %                 else
                %                     pRapAbAg = exp(-obj.EAgMem)/Z;
                %                 end
                %                 for i=1:BitOn
                %                     if(rand<pRapAbAg)
                %                         TotAg = TotAg+1;
                %                     end
                %                 end
                
                if(rand<Pcapture)
                    TotAg = TotAg+1;
                end
                
                %                 for i=1:BitOn
                %                     switch obj.BindFlag
                %                         case 0
                %                             if(rand<Pcapture)
                %                                 TotAg = TotAg+1;
                %                             end
                %                         case 1
                %                             if(rand<Pcapture)
                %                                 TotAg = TotAg+1;
                %                             end
                %                         case 3
                %                             if(rand<Pcapture)
                %                                 TotAg = TotAg+1;
                %                             end
                %                     end
                %                 end
            end
        end
        
        %%%
        %%%%%%
        %%%%%%
        %%%%%%
        
        function TotAg = CAgCaptureDensity2Arms2Ags(obj);
            
            kb = 1.38064852e-23;
            %             F = 20e-12;
            %             F = 1e-12;
            F = obj.F;
            
            T=273;
            xb = 1e-9;
            xb*F/(kb*T);
            
            
            
            
            lambda = exp(-obj.EAgMem)*exp(xb*F/(kb*T));
            r = exp(-obj.EAbAg)*exp(xb*F/(kb*T));
            k = exp(-obj.EAbAg)*exp(0.5*xb*F/(kb*T));
            xi =  exp(-obj.EAgMem)*exp(0.5*xb*F/(kb*T));
            
            q = (obj.NAg-1)*obj.q;
            
            %             switch obj.BindFlag
            %                 case 1
            %                     q = (obj.NAg-1)*obj.q;
            %                 case 2
            %                     q = (obj.NAg-1)*obj.q;
            %                 case 3
            %                     q = obj.q;
            %             end
            
            nu = 0;
            
            %           At t=0, either bond arm are bound/one of the arm is bound/ none
            %           of the arms is bound
            Cap15_2 = lambda.*xi.*(q+r+xi).*(k.^2.*(r+xi)+lambda.*xi.*(q+r+xi)+k.*(r.*(lambda+nu+xi)+xi.*(q+lambda+nu+xi))).^(-1);
            
            %             lambda.^2.*(q+r+lambda).*(k.^2.*(r+lambda)+lambda.^2.*(q+r+lambda)+k.*(q.*lambda+(r+lambda).*(2.*lambda+nu) )).^(-1);
            
            Cap15_1 = k.*xi.*(k+q+r+lambda+nu+xi).*(k.^2.*(r+xi)+lambda.*xi.*(q+r+xi)+k.*(r.*(lambda+nu+xi)+xi.*(q+lambda+nu+xi))).^(-1);
            
            %             k.*lambda.*(k+q+r+2.*lambda+nu).*(k.^2.*(r+lambda)+ lambda.^2.*(q+r+lambda)+k.*(q.*lambda+(r+lambda).*( 2.*lambda+nu))).^(-1);
            
            Cap14_13_2 = (lambda*(k*(r + xi) + xi*(q + r + xi)))/(k^2*(r + xi) + lambda*xi*(q + r + xi) +  k*(r*(lambda + nu + xi) + xi*(q + lambda + nu + xi)));
            
            %             lambda.*(k.*(r+lambda)+lambda.*(q+r+lambda)).*(k.^2.*(r+lambda)+lambda.^2.*(q+r+lambda)+k.*(q.*lambda+(r+lambda).*(2.*lambda+nu))).^(-1);
            %
            Cap14_13_1 = (k*(k*(r + xi) + xi*(q + r + nu + xi)))/(k^2*(r + xi) + lambda*xi*(q + r + xi) +  k*(r*(lambda + nu + xi) + xi*(q + lambda + nu + xi)));
            
            
            %             k.*(k.*(r+lambda)+lambda.*(q+r+lambda+nu)).*(k.^2.*(r+lambda)+lambda.^2.*(q+r+lambda)+k.*(q.*lambda+(r+lambda).*(2.*lambda+nu))).^(-1);
            
            Cap11_7_2 = (q*lambda*xi)/(k^2*(r + xi) + lambda*xi*(q + r + xi) + k*(r*(lambda + nu + xi) + xi*(q + lambda + nu + xi)));
            
            %             q.*lambda.^2.*(k.^2.*(r+lambda)+lambda.^2.*(q+r+lambda)+k.*(q.*lambda+(r+lambda).*(2.*lambda+nu))).^(-1);
            
            Cap11_7_1 = xi.*(k.^2+lambda.*xi+k.*(q+lambda+nu+xi)).*(k.^2.*(r+xi)+lambda.*xi.*(q+r+xi)+k.*(r.*(lambda+nu+xi)+xi.*(q+lambda+nu+xi))).^(-1);
            
            %             lambda.*(k.^2+lambda.^2+k.*(q+2.*lambda+nu)).*(k.^2.*(r+lambda)+lambda.^2.*(q+r+lambda)+k.*(q.*lambda+(r+lambda).*(2.*lambda+nu))).^(-1);
            
            NBCR = obj.NBCR;
            TotAg = 0;
            P = cumsum(obj.Parr);
            for i=1:NBCR
                ii = min(find(P-rand>0));
                idx = strfind(dec2bin(ii-1,obj.BondsNum),'1');
                BitOn = length(idx);
                switch BitOn
                    case 0
                        %                         E=0;
                        Pcapture = 0;
                    case 1 % starting from 11 or 7
                        Q = cumsum([Cap11_7_1 Cap11_7_2 (1-Cap11_7_1-Cap11_7_2)]);
                        ii = min(find(Q-rand>0));
                        if(ii==1)
                            TotAg = TotAg+1;
                        elseif(ii==2)
                            TotAg = TotAg+2;
                        end
                    case 2 % starting from 15
                        Q = cumsum([Cap15_1 Cap15_2 (1-Cap15_1-Cap15_2)]);
                        ii = min(find(Q-rand>0));
                        if(ii==1)
                            TotAg = TotAg+1;
                        elseif(ii==2)
                            TotAg = TotAg+2;
                        end
                end
                
            end
        end
        
        %%%
        %%%
        %%%
        
        function TotAg = CAgCaptureDensityCluster(obj);
            
            ClusterSize = obj.AbAgClusterProp.ClusterSize;
            NClusters = round(obj.NBCR/ClusterSize);
            TotAg = 0;
            
            %%%%%%
            %%%%%%
            %%%%%%
            if(obj.CapSimFlag==0)
                if(obj.EqFlag==0)
                    AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                    P = cumsum(obj.Parr);
                    for i=1:NClusters
                        ii = min(find(P-rand>0));
                        switch ii
                            case 1
                                %                         E=0;
                                Pcapture = 0;
                            otherwise
                                Q = cumsum(AgCapProStartingState,2);
                                ii2 = min(find(Q(ii-1,:)-rand>0)) - 1; % Amount of Ag added
                                TotAg = TotAg + ii2;
                        end
                    end
                else
                    
                    AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                    P = cumsum(obj.Parr);
                    N = length(P);
                    for i=1:NClusters
                        ii = min(find(P-rand>0));
                        switch ii
                            case N
                                %                         E=0;
                                Pcapture = 0;
                            otherwise
                                Q = cumsum(AgCapProStartingState,2);
                                ii2 = min(find(Q(ii,:)-rand>0)) - 1; % Amount of Ag added
                                TotAg = TotAg + ii2;
                        end
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%
            else
                
                P = cumsum(obj.AbAgClusterProp.Parr);
                
                countArr = [];
                EndStateNames = {zeros(1,ClusterSize)};
                EndStateNamesHits = [0];
                AgCapDist = zeros(1,2*ClusterSize+1);
                CaptureTimes = 0;
                for j=1:NClusters
                    ii = min(find(P-rand>0));
                    if(ii>1)
                        CaptureTimes = CaptureTimes+1;
                        InitialState = obj.AbAgClusterProp.StartingStateNames{ii-1};
                        [EndStateName EndStateNum ACstates_AgCap count]  = obj.AbAgClusterProp.SimulateClusterAbs(InitialState);
                        I = ismember(cell2mat(EndStateNames),EndStateName,'rows');
                        idxAlready = find(I);
                        if(length(idxAlready)==0)
                            EndSize = length(EndStateNames);
                            EndStateNames{EndSize+1,1} = EndStateName;
                            EndStateNamesHits(EndSize+1,1) = 1;
                        else
                            EndStateNamesHits(idxAlready,1) = EndStateNamesHits(idxAlready,1)+1;
                        end
                        AgCapDist(ACstates_AgCap+1) = AgCapDist(ACstates_AgCap+1) + 1;
                        countArr(j) = count;
                        TotAg = TotAg + ACstates_AgCap;
                    end
                end
                AgCapDist/CaptureTimes;
                countArr(countArr>1);
            end
            
        end
        
        %%%
        %%%
        %%%
        
        %%%
        %%%
        %%%
        
        function [TotAg1 TotAg2] = CAgCaptureDensityCluster2Ag(obj);
            
            global ICArr
            
            ClusterSize = obj.AbAgClusterProp.ClusterSize;
            NClusters = round(obj.NBCR/ClusterSize);
            TotAg1 = 0;
            TotAg2 = 0;
            
            %%%%%%
            %%%%%%
            %%%%%%
            if(obj.CapSimFlag==0)
                if(obj.EqFlag==0)
                    for n_bcr = 1:obj.aff.NBCR
                        AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingStateArr{n_bcr};
                        %                     AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                        %                     P = cumsum(obj.Parr);
                        P = cumsum(obj.AbAgClusterProp.ParrNBCR{n_bcr});
                        %                     for i=1:NClusters
                        randomnumber1 = rand;
                        ii = min(find(P-randomnumber1>0));
                        if( (isnumeric(ii)==0) || (length(ii)>1) )
                            obj.AbAgClusterProp.AgCapProStartingStateArr{n_bcr}
                            obj.AbAgClusterProp.ParrNBCR{n_bcr}
                            P
                            randomnumber1
                            ii
                        end
                        switch ii
                            case 1
                                %                         E=0;
                                Pcapture = 0;
                            otherwise % ii chooses amongst the different states in obj.AbAgClusterProp.StartingStateNames
                                Q = cumsum(AgCapProStartingState,2);
                                % now choose an amount of captured Ag from obj.AbAgClusterProp.AgCapState
                                randnumber = rand;
                                ii2 = min(find(round(Q(ii-1,:),4)-randnumber>=0));
                                %                                 ii2 = min(find(Q(ii-1,:)-randnumber>0));
                                %                                 ii2 = min(find(Q(ii-1,:)-rand>0));
                                CapturedAgs = obj.AbAgClusterProp.AgCapState(ii2,:);
                                
                                
                                
                                %                                 if( CapturedAgs(1)>0 | CapturedAgs(2) > 0) % remove the captured Ags from the immune complex (IC)
                                %                                     IC_encounterd = obj.AbAgClusterProp.ICnum;
                                %                                     ICremoveFlag = 1;
                                %                                 else
                                %                                     ICremoveFlag = 0;
                                %
                                %                                 end
                                if(obj.aff.ExtcFlag==1) % extract the whole immune complex
                                    if(size(CapturedAgs,1)==0)
                                        n_bcr
                                        AgCapProStartingState
                                        obj.AbAgClusterProp.ParrNBCR{n_bcr}
                                        Q
                                        randnumber
                                        ii2
                                        obj.AbAgClusterProp.AgCapState
                                        CapturedAgs
                                        continue;
                                    end
                                    if( CapturedAgs(1)>0 | CapturedAgs(2) > 0) % remove the captured Ags from the immune complex (IC)
                                        IC_encounterd = obj.AbAgClusterProp.ICnum;
                                        TotAg1 = TotAg1 + ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}(1);
                                        TotAg2 = TotAg2 + ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}(2);
                                        ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr} = single([0 0]);
                                    end
                                elseif(obj.aff.ExtcFlag==2)
                                    %                                 Removing the captured Ag from the IC
                                    TotAg1 = TotAg1 + CapturedAgs(1);
                                    TotAg2 = TotAg2 + CapturedAgs(2);
                                    if( CapturedAgs(1)>0 | CapturedAgs(2) > 0) % remove the captured Ags from the immune complex (IC)
                                        IC_encounterd = obj.AbAgClusterProp.ICnum;
                                        temp = ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr} - CapturedAgs;
                                        if(temp(1) <0 | temp(2) <0)
                                            temp
                                        end
                                        ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr} = ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr} - CapturedAgs;
                                        if(ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}(1)<0 | ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}(2)<0)
                                            ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}
                                        end
                                    end
                                end
                                %                                 ii2 = min(find(Q(ii-1,:)-rand>0)) - 1; % Amount of Ag added
                                %                                 TotAg = TotAg + ii2;
                        end
                    end
                    
                    
                    if(obj.aff.ExtcFlag==1)
                        IC_encounterd = obj.AbAgClusterProp.ICnum;
                        emptyflag = obj.isemptyIC(IC_encounterd);
                        %                         if(emptyflag)
                        %                             ICArr{IC_encounterd(1)}(IC_encounterd(2)) = [];
                        %                         end
                    end
                    
                    %                     if(ICremoveFlag)
                    %                         idx_IC = [1:length(ICArr{IC_encounterd(1)})];
                    %                         idx2remove = IC_encounterd(2);
                    %                         ICArr_temp = setdiff(idx_IC,idx2remove);
                    %
                    %                         ICArr{IC_encounterd(1)} = ICArr{IC_encounterd(1)}{ICArr_temp};
                    %                     end
                    
                else
                    
                    AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                    P = cumsum(obj.Parr);
                    N = length(P);
                    for i=1:NClusters
                        ii = min(find(P-rand>0));
                        switch ii
                            case N
                                %                         E=0;
                                Pcapture = 0;
                            otherwise
                                Q = cumsum(AgCapProStartingState,2);
                                ii2 = min(find(Q(ii,:)-rand>0)) - 1; % Amount of Ag added
                                TotAg = TotAg + ii2;
                        end
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%
            else
                
                P = cumsum(obj.AbAgClusterProp.Parr);
                
                countArr = [];
                EndStateNames = {zeros(1,ClusterSize)};
                EndStateNamesHits = [0];
                AgCapDist = zeros(1,2*ClusterSize+1);
                CaptureTimes = 0;
                for j=1:NClusters
                    ii = min(find(P-rand>0));
                    if(ii>1)
                        CaptureTimes = CaptureTimes+1;
                        InitialState = obj.AbAgClusterProp.StartingStateNames{ii-1};
                        [EndStateName EndStateNum ACstates_AgCap count]  = obj.AbAgClusterProp.SimulateClusterAbs(InitialState);
                        I = ismember(cell2mat(EndStateNames),EndStateName,'rows');
                        idxAlready = find(I);
                        if(length(idxAlready)==0)
                            EndSize = length(EndStateNames);
                            EndStateNames{EndSize+1,1} = EndStateName;
                            EndStateNamesHits(EndSize+1,1) = 1;
                        else
                            EndStateNamesHits(idxAlready,1) = EndStateNamesHits(idxAlready,1)+1;
                        end
                        AgCapDist(ACstates_AgCap+1) = AgCapDist(ACstates_AgCap+1) + 1;
                        countArr(j) = count;
                        TotAg = TotAg + ACstates_AgCap;
                    end
                end
                AgCapDist/CaptureTimes;
                countArr(countArr>1);
            end
            
        end
        
        %%%
        %%%
        %%%
        
        
        function TotAg = CAgCaptureDensityClusterDistribution(obj);
            
            ClusterSize = obj.AbAgClusterProp.ClusterSize;
            NClusters = round(obj.NBCR/ClusterSize);
            TotAg = 0;
            
            %%%%%%
            %%%%%%
            %%%%%%
            if(obj.CapSimFlag==0)
                if(obj.EqFlag==0)
                    AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                    P = cumsum(obj.Parr);
                    for i=1:NClusters
                        
                        %%% randomize amount of Ag
                        
                        obj.AgDistribution(obj.MAg);
                        if(obj.CapSimFlag)
                            obj.AbAgClusterProp = SimulateClusterAbsClass(obj.NAg,obj.q,obj.EAgMem,obj.EAbAg,obj.tau,obj.F,obj.AbAgClusterProp.ClusterSize);
                        else
                            obj.AbAgClusterProp = AbAgCluster(aff,obj.NBCR,obj.NAg,obj.q,obj.EAbAg,obj.tau,obj.F,obj.AbAgClusterProp.ClusterSize,obj.EqFlag);
                        end
                        obj.Parr = obj.AbAgClusterProp.Parr;
                        
                        AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                        P = cumsum(obj.Parr);
                        
                        %%%%%%
                        
                        
                        
                        
                        
                        
                        ii = min(find(P-rand>0));
                        switch ii
                            case 1
                                %                         E=0;
                                Pcapture = 0;
                            otherwise
                                Q = cumsum(AgCapProStartingState,2);
                                ii2 = min(find(Q(ii-1,:)-rand>0)) - 1; % Amount of Ag added
                                TotAg = TotAg + ii2;
                        end
                    end
                else
                    
                    AgCapProStartingState = obj.AbAgClusterProp.AgCapProStartingState;
                    P = cumsum(obj.Parr);
                    N = length(P);
                    for i=1:NClusters
                        ii = min(find(P-rand>0));
                        switch ii
                            case N
                                %                         E=0;
                                Pcapture = 0;
                            otherwise
                                Q = cumsum(AgCapProStartingState,2);
                                ii2 = min(find(Q(ii,:)-rand>0)) - 1; % Amount of Ag added
                                TotAg = TotAg + ii2;
                        end
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%
            else
                
                P = cumsum(obj.AbAgClusterProp.Parr);
                
                countArr = [];
                EndStateNames = {zeros(1,ClusterSize)};
                EndStateNamesHits = [0];
                AgCapDist = zeros(1,2*ClusterSize+1);
                CaptureTimes = 0;
                for j=1:NClusters
                    ii = min(find(P-rand>0));
                    if(ii>1)
                        CaptureTimes = CaptureTimes+1;
                        InitialState = obj.AbAgClusterProp.StartingStateNames{ii-1};
                        [EndStateName EndStateNum ACstates_AgCap count]  = obj.AbAgClusterProp.SimulateClusterAbs(InitialState);
                        I = ismember(cell2mat(EndStateNames),EndStateName,'rows');
                        idxAlready = find(I);
                        if(length(idxAlready)==0)
                            EndSize = length(EndStateNames);
                            EndStateNames{EndSize+1,1} = EndStateName;
                            EndStateNamesHits(EndSize+1,1) = 1;
                        else
                            EndStateNamesHits(idxAlready,1) = EndStateNamesHits(idxAlready,1)+1;
                        end
                        AgCapDist(ACstates_AgCap+1) = AgCapDist(ACstates_AgCap+1) + 1;
                        countArr(j) = count;
                        TotAg = TotAg + ACstates_AgCap;
                    end
                end
                AgCapDist/CaptureTimes;
                countArr(countArr>1);
            end
            
        end
        
        %%%
        %%%
        %%%
        
        
        %%%
        %%%
        %%%
        function AgDistribution(obj,MAg)
            obj.MAg = MAg;
            R = 100;
            A = obj.A;
            p = (A/R)^2;
            %             p=0.1;
            obj.NAg = binornd(MAg,p);
        end
        
        
        
        function AgDistributionNanoParticle(obj,MAg)
            % We assume that the BCR has access to all antigens
%             obj.MAg = MAg;
%             R = 100;
%             A = obj.A;
%             p = (A/R)^2;
%             %             p=0.1;
%             obj.NAg = binornd(MAg,p);
            obj.NAg = MAg;
        end          
        
        function emptyflag = isemptyIC(obj,IC_encounterd)
            global ICArr
            emptyflag = 1;
            for n_bcr=1:length(ICArr{IC_encounterd(1)}{IC_encounterd(2)})
                for j=1:length(length(ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}))
                    if(ICArr{IC_encounterd(1)}{IC_encounterd(2)}{n_bcr}(j))
                        emptyflag = 0;
                    end
                end
            end
        end
      
        
    end
    
end

