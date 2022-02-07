classdef Bcell < handle
    properties
        BCR_s;
        AffinityClass;
        lambda;
        mu;
        Diff;
        dt;
        x;
        X;
        RepTime;
        DeathTime;
        T2rep;
        Age;
        N;
        ParentPos;
        Sons;
        SelfNum;
        ParentNum;
        tauNonlinear=0.5;
        timeln;
        Flagln;
        lambda0;
        Dev0;
        DevLeft;
        alpha;
        Jackpotflag;
        Nmax;
        lambdaBase; %% factor for affinitny to birthrate
        lambdaBaseExp; %% Lambda for B cell expansion outside of the GC       
        AgCap;
        BondsNum;
        NBCR;
        SecBind;
        paramFlag;
        AgCapm1;
        AgCapm2;
        AgCapmArr;
        Mature;
        lambdaPersonal;
        %         Sons;
        %         TotAgents;
        %         nm;
        %         len;
    end
    methods
        %         function obj = Bcell(lambda,mu,Diff,xinit,dt,ParentNum,Jackpotflag,JPThreshold,Nmax,lambdaBase,BondsNum,NBCR);
        function obj = Bcell(BCR_s);
            if(nargin > 0)
                %                 BCR_s = struct('l',lambda,'mu',mu,'D',Diff,'n',n-1,'dt',dt,'ParentNum',0,'Jackpotflag',Jackpotflag,'JPThreshold',JPThreshold,'Nmax',Nmax,'lambdaBase',lambda0,'BondsNum',BondsNum,'NBCR',NBCR);
                %                 obj.x = xinit;
                
                obj.BCR_s = BCR_s;
                obj.X = BCR_s.n;
                obj.N = length(obj.X);
                obj.lambda = BCR_s.l;
                obj.mu = BCR_s.mu;
                obj.Diff = BCR_s.D;
                obj.dt = BCR_s.dt;
                obj.T2rep = 0;
                obj.Age = 0;
                obj.ParentPos = BCR_s.n;
                obj.SelfNum =  rand;
%                 obj.ParentNum = BCR_s.ParentNum;
                
                obj.Mature = 0;
                
                obj.Flagln=0;
                obj.lambda0=BCR_s.JPThreshold;
                obj.Dev0=4;
                obj.alpha=2;
                obj.Jackpotflag=BCR_s.Jackpotflag;
                obj.Nmax = BCR_s.Nmax;
                obj.lambdaBase = BCR_s.lambdaBase;
                obj.lambdaBaseExp = BCR_s.lambdaBaseExp;
                %                 obj.BondsNum = BCR_s.BondsNum;
                %                 obj.NBCR = BCR_s.NBCR;
                %                 obj.SecBind = BCR_s.SecBind;
                %                 obj.paramFlag = BCR_s.paramFlag;
                %                 a = Affinity(obj.lambda,obj.paramFlag,obj.BondsNum,obj.NBCR,obj.SecBind);
                %                 a = Affinity(obj.lambda,BCR_s.aff);
                obj.AffinityClass = Affinity(obj.lambda,BCR_s.aff);
                %                 obj.AgCap = a.CAgCapture();
                %                 obj.AgCap = a.CAgCaptureDensity();
                
                obj.AgCap = obj.AffinityClass.CAgCaptureFunc();
                
                %                 obj.AgCaptureInit();
                
                %                 obj.AgCapm1=1e3;
                %                 obj.AgCapm2=1e3;
                
                %                 obj.BCR_s = BCR_s;
                %                 obj.X = xinit;
                %                 obj.N = length(obj.X);
                %                 obj.lambda = lambda;
                %                 obj.mu = mu;
                %                 obj.Diff = Diff;
                %                 obj.dt = dt;
                %                 obj.T2rep = 0;
                %                 obj.Age = 0;
                %                 obj.ParentPos = xinit;
                %                 obj.SelfNum =  rand;
                %                 obj.ParentNum = ParentNum;
                %
                %                 obj.Flagln=0;
                %                 obj.lambda0=JPThreshold;
                %                 obj.Dev0=4;
                %                 obj.alpha=2;
                %                 obj.Jackpotflag=Jackpotflag;
                %                 obj.Nmax = Nmax;
                %                 obj.lambdaBase = lambdaBase;
                %                 obj.BondsNum = BondsNum;
                %                 obj.NBCR = NBCR;
                % %                 obj.AgCap = 10;
                %                 a = Affinity(obj.lambda,obj.BondsNum,obj.NBCR);
                %                 obj.AgCap = a.CAgCapture();
            else
            end
        end
        function move(obj);
            A = obj.dt;
            noise_dev = sqrt(2*A);
            f = normrnd(0,noise_dev);
            dx = f;
            obj.x = obj.x + dx;
        end
        function MoveLattice(obj);
            N = obj.N;
            p = obj.dt*obj.D;
            dX = zeros(N,1);
            
            if(rand <= p)
                dX = zeros(N,1);
                k = ceil(N*rand);
                if(k==0)
                    k
                end
                if(rand <= 0.5)
                    dX(k) = dX(k)+1;
                else
                    dX(k) = dX(k)-1;
                end
                obj.X = obj.X + dX;
            end
        end
        function MoveLatticeUp(obj);
            N = obj.N;
            p = obj.dt*obj.D;
            dX = zeros(N,1);
            
            if(rand <= p)
                dX = zeros(N,1);
                k = ceil(N*rand);
                if(k==0)
                    k
                end
                dX(k) = dX(k)+1;
                obj.X = obj.X + dX;
            end
        end
        function MoveRandJump(obj);
            N = obj.N;
            p = obj.dt*obj.D;
            dX = zeros(N,1);
            
            if(rand <= p)
                dX = zeros(N,1);
                k = ceil(N*rand);
                if(k==0)
                    k
                end
                dX(k) = dX(k)+10*rand;
                obj.X = obj.X + dX;
            end
        end
        function MoveDiffusion(obj);
            N = obj.N;
            p = obj.dt*obj.D;
            dX = zeros(N,1);
            
            if(rand <= p)
                dX = zeros(N,1);
                k = ceil(N*rand);
                if(k==0)
                    k
                end
                
                A = obj.dt;
                noise_dev = sqrt(2*A);
                dX(k) = normrnd(0,noise_dev);
                
                obj.X = obj.X + dX;
            end
        end
        function Kid = RepCheck(obj);
            obj.T2rep = obj.T2rep + obj.dt;
            if(obj.T2rep>=obj.RepTime)
                Kid = BC(obj.lambda,obj.mu,obj.x,obj.dt);
                obj.T2rep = 0;
                obj.RepTime=random('Exponential',1/obj.lambda);
            else
                Kid = [];
            end
        end
        function Dflag = DeathCheck(obj);
            obj.Age = obj.Age + obj.dt;
            if(obj.Age>=obj.DeathTime)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        
        function Kid = Birth(obj)
            if(obj.Jackpotflag)
                Kid = obj.BirthJackpotKidMut();
            else
                Kid = obj.BirthKidMut();
            end
        end
        
        function Kid = BirthNorm(obj,lambdaMeanAll,CTMeanAll)
            if(obj.Jackpotflag)
                Kid = obj.BirthJackpotKidMut();
            else
                Kid = obj.BirthKidMutNorm(lambdaMeanAll,CTMeanAll);
            end
        end
        
        
        %%%%%%%%
        
        function Kid = BirthNormMemoryOutsideGC(obj,lambdaMeanAll,CTMeanAll)
%             AgThExp = 100;
            AgThExp = obj.BCR_s.aff.AgThExp;

             CT = obj.BCR_s.aff.CT; 
             %CT = 000; % Does it depend on TfhCs? Sould not not
	    CTMeanAll=CT;
            if(obj.AgCap>=AgThExp)
                if(lambdaMeanAll)
                    lambda = obj.BCR_s.b*obj.BCR_s.lambdaBaseExp*( (CT + obj.AgCap)/(CTMeanAll + lambdaMeanAll));
                else
                    lambda = obj.BCR_s.b*obj.BCR_s.lambdaBaseExp;
                end
            else
                lambda = 0;
            end
            obj.lambdaPersonal = lambda;
            if(isnan(lambda)) lambda=0; end
            if(length(lambda)>1)
                Kid = [];
                return;
            end
            
            p = obj.dt*lambda;
            if(rand <= p)

                BCR_s_Kid = obj.BCR_s;
                FlagMutOk = 0;
                
                rcv_kid = BCR_s_Kid.aff.rcv;

                
                EAbAg_new = obj.BCR_s.aff.EAbAg;

                if( EAbAg_new > -1)
                    BCR_s_Kid.aff.EAbAg = EAbAg_new;
                    FlagMutOk = 1;
                end

                global mutrate
                if isempty(mutrate) %makes it so mutrate is no longer needed for other initializaiton files
                else
                if rand <mutrate %introduce mutaiton to B cells.
                            MutStd = sqrt(2*obj.Diff);
                            CorrVC = obj.BCR_s.aff.CorrVC; % the correlation between the effect of a mutation in the concerved site on the binding energy to the variable domain
                            CorrVV = obj.BCR_s.aff.CorrVV; % the correlation between the effect of a mutation in the variable site on the binding energy to the concerved domain
                            MutVC = obj.BCR_s.aff.MutVC; % the independent change in enrgy following a mutation in the concerved site on the binding energy to the variable domain
                            MutVV = obj.BCR_s.aff.MutVV; %

                            Ec = obj.BCR_s.aff.EAbAg(1);       Ev1 = obj.BCR_s.aff.EAbAg(2);
                            Ev2 = obj.BCR_s.aff.EAbAg(3);       Ev3 = obj.BCR_s.aff.EAbAg(4);
      
             pvar=(copularnd('t',CorrVV,1,1));  %replace mutation scheme with copular to 

  %dE=zeros(length(pvar),1);
%dE(pvar<=0.0793848)=log(pvar(pvar<=0.0793848)/0.0793848)/1.43256;
%dE(pvar >0.0793848)=log(pvar((pvar>0.0793848))/0.920617)/-0.329412;
%dE=-dE;
             mu=1;                          %generate correlated uniform distributions
            sigma=.5;     offset=4;
            dE = MutStd*(logninv(pvar,mu,sigma)-offset);  %take value of correlations then do inverse function to find value
            Ev1=Ev1+dE(1);     Ev2=Ev2+dE(2);    Ev3=Ev3+dE(3);

                            Ec = max(Ec,-1);          Ev1 = max(Ev1,-2);
                            Ev2 = max(Ev2,-2);   Ev3 = max(Ev3, -2);
                            BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                            BCR_s_Kid.aff.EAbAg(1) = Ec;   BCR_s_Kid.aff.EAbAg(2) = Ev1;
                            BCR_s_Kid.aff.EAbAg(3) = Ev2;   BCR_s_Kid.aff.EAbAg(4) = Ev3;
                else
                end
                end           
                if(FlagMutOk)
                    Kid = Bcell(BCR_s_Kid);
                    
                    Kid.Mature = 1; % This is a proliferation of memory B cells. The kid is mature (memory)
                    
                    obj.Sons = [obj.Sons ; Kid.SelfNum];
                else
                    Kid = [];
                end
            else
                Kid = [];
            end
      
        end
        
        %%%%%%%%
        
        function Kid = BirthT(obj)
            lambda = obj.lambda;
            %             x = obj.x;
            
            p = obj.dt*lambda;
            
            if(rand <= p)
                Kid = BDD_moduleT(obj.lambda,obj.mu,obj.D,obj.X,obj.dt,obj.SelfNum);
                obj.Sons = [obj.Sons ; Kid.SelfNum];
            else
                Kid = [];
            end
        end
        
        %%%%%%%%%%%%%%
        
        function Kid = BirthTX(obj,lambdaX)
            
            %             X = obj.X;
            %             Y = int8(X'+1);
            %             cellJ = num2cell(Y);
            lambda = obj.lambdaBase;
            
            p = obj.dt*lambda;
            if(rand <= p)
                BCR_s_Kid = obj.BCR_s;
                %                 Kid = Bcell(obj.lambda,obj.mu,obj.Diff,obj.X,obj.dt,obj.SelfNum,obj.Jackpotflag,obj.lambda0,obj.Nmax,obj.lambdaBase,obj.BondsNum,obj.NBCR);
                Kid = Bcell(BCR_s_Kid);
                
                obj.Sons = [obj.Sons ; Kid.SelfNum];
            else
                Kid = [];
            end
        end
        
        %%%%%%%%%%%%%
        
        function Kid = BirthTXNonLinear(obj,lambdaX)
            
            X = obj.X;
            
            Y = int8(X'+1);
            cellJ = num2cell(Y);
            
            lambda = obj.lambda;
            
            if(length(lambda)>1)
                Kid = [];
                return;
            end
            if(obj.Flagln)
                lambda = lambda*(1+exp(-obj.timeln/obj.tauNonlinear));
                obj.timeln = obj.timeln+0.01;
            end
            if(  (lambda>obj.lambda0) && (obj.Flagln==0))
                obj.timeln=0;
                obj.Flagln=1;
                lambda = lambda*(1+exp(-obj.timeln/obj.tauNonlinear));
            end
            
            if(obj.timeln>1)
                obj.Flagln = 0;
            end
            
            p = obj.dt*lambda;
            
            if(rand <= p)
                %                 Kid = Bcell(obj.lambda,obj.mu,obj.Diff,obj.X,obj.dt,obj.SelfNum,obj.Jackpotflag,obj.lambda0,obj.Nmax,obj.lambdaBase,obj.BondsNum,obj.NBCR);
                BCR_s_Kid = obj.BCR_s;
                Kid = Bcell(BCR_s_Kid);
                obj.Sons = [obj.Sons ; Kid.SelfNum];
            else
                Kid = [];
            end
        end
        
        %%%%%%%%%%%%%
        
        
        function Kid = BirthKidMut(obj)
            
            lambda = obj.lambda;
            
            if(length(lambda)>1)
                Kid = [];
                return;
            end
            
            p = obj.dt*lambda;
            if(rand <= p)
                %                 MutStd = sqrt(2*obj.Diff*obj.dt);
                MutStd = sqrt(2*obj.Diff);
                lambdaSon = abs(obj.lambda + normrnd(0,MutStd));
                %                 Kid = Bcell(lambdaSon,obj.mu,obj.Diff,obj.X,obj.dt,obj.SelfNum,obj.Jackpotflag,obj.lambda0,obj.Nmax,obj.lambdaBase,obj.BondsNum,obj.NBCR);
                BCR_s_Kid = obj.BCR_s;
                BCR_s_Kid.l = lambdaSon;
                Kid = Bcell(BCR_s_Kid);
                obj.Sons = [obj.Sons ; Kid.SelfNum];
            else
                Kid = [];
            end
        end
        
        %%%%%%%%%%%%%
        
        function Kid = BirthKidMutNorm(obj,lambdaMeanAll,CTMeanAll)
            CT = obj.BCR_s.aff.CT;
            %             lambda = obj.lambdaBase*obj.AgCap/lambdaMeanAll;
            %             lambda = obj.lambdaBase*( (CT + obj.AgCap)/(CT + lambdaMeanAll));
            if(obj.AgCap>0)
%                 lambda = obj.lambdaBase*( (CT + obj.AgCap)/(CT + lambdaMeanAll));
%                 lambda = obj.BCR_s.b*obj.lambdaBase*( (CT + obj.AgCap)/(CT + lambdaMeanAll));
            if( (lambdaMeanAll==0) && (CTMeanAll==0))
                lambda = obj.BCR_s.b*obj.lambdaBase;
            else
global expon

                lambda = obj.BCR_s.b*obj.lambdaBase*( (CT^expon + obj.AgCap^expon)/(CTMeanAll^expon + lambdaMeanAll));
            end
                
            else
                lambda = 0;;
            end
            obj.lambdaPersonal = lambda;
            if(isnan(lambda)) lambda=0; end
            if(length(lambda)>1)
                Kid = [];
                return;
            end
            
            p = obj.dt*lambda;
            if(rand <= p)
                %                 MutStd = sqrt(2*obj.Diff*obj.dt);
                MutStd = sqrt(2*obj.Diff);
                BCR_s_Kid = obj.BCR_s;
                FlagMutOk = 0;
                
                rcv_kid = BCR_s_Kid.aff.rcv;
                %                 if(rand<0.1)
                %                     if(rand<0.5)
                %                         rcv_kid(1) = abs(rcv_kid(1) + normrnd(0,MutStd));
                %                     else
                %                         rcv_kid(2) = abs(rcv_kid(2) + normrnd(0,MutStd));
                %                     end
                %                     rcv_kid = rcv_kid/sum(rcv_kid);
                %                 end
                %                 BCR_s_Kid.aff.rcv = rcv_kid;
                
                switch obj.BCR_s.aff.paramFlag
                    case 1
                        if(obj.BCR_s.aff.BindFlag==9)
                            if(rand<obj.BCR_s.aff.wcv(1)) %mutation in the conserved part
                                EAbAg_cons_new = obj.BCR_s.aff.EAbAg(1) + normrnd(0,MutStd);
                                if(EAbAg_cons_new<0) EAbAg_cons_new=0; end
                                if( EAbAg_cons_new >= 0)
                                    BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                                    %                                     BCR_s_Kid.aff.rcv = obj.ModifyRCV(obj.BCR_s.aff.EAbAg(1),EAbAg_cons_new,rcv_kid,1);
                                    BCR_s_Kid.aff.EAbAg(1) = EAbAg_cons_new;
                                    FlagMutOk = 1;
                                end
                            else %mutation in the variable part ; We take them to have a correlation of -2D such that increase for antigen
                                % a results in a decrease for b
                                dE = normrnd(0,MutStd);
                                EAbAg_var1_new = obj.BCR_s.aff.EAbAg(2) + dE;
                                EAbAg_var2_new = obj.BCR_s.aff.EAbAg(3) - dE;
                                %                                 EAbAg_var1_new = obj.BCR_s.aff.EAbAg(2) + normrnd(0,MutStd);
                                %                                 EAbAg_var2_new = obj.BCR_s.aff.EAbAg(3) + normrnd(0,MutStd);
                                if(EAbAg_var1_new<0) EAbAg_var1_new=0; end
                                if(EAbAg_var2_new<0) EAbAg_var2_new=0; end
                                if( ( EAbAg_var1_new >= 0) & ( EAbAg_var2_new >= 0) )
                                    BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                                    BCR_s_Kid.aff.EAbAg(2) = EAbAg_var1_new;
                                    BCR_s_Kid.aff.EAbAg(3) = EAbAg_var2_new;
                                    FlagMutOk = 1;
                                end
                            end
                        else
                            EAbAg_new = obj.BCR_s.aff.EAbAg + normrnd(0,MutStd);
                            if( EAbAg_new > 0)
                                BCR_s_Kid.aff.EAbAg = EAbAg_new;
                                FlagMutOk = 1;
                            end
                        end
                        %                         BCR_s_Kid.aff.EAbAg = abs(obj.BCR_s.aff.EAbAg + normrnd(0,MutStd));
                    case 2
                        
                        if(obj.BCR_s.aff.BindFlag==9)
                            if(rand<obj.BCR_s.aff.wcv(1)) %mutation in the conserved part
                                q1_cons_new = obj.BCR_s.aff.q1(1) + normrnd(0,MutStd);
                                if(q1_cons_new<0) q1_cons_new=0; end
                                if( q1_cons_new >= 0)
                                    BCR_s_Kid.aff.q1 = obj.BCR_s.aff.q1;
                                    BCR_s_Kid.aff.q1(1) = q1_cons_new;
                                    FlagMutOk = 1;
                                end
                            else %mutation in the variable part ; We take them to have a correlation of -2D such that increase for antigen
                                % a results in a decrease for b
                                dq = normrnd(0,MutStd);
                                q1_var1_new = obj.BCR_s.aff.q1(2) + dq;
                                q1_var2_new = obj.BCR_s.aff.q1(3) - dq;
                                %                                 q1_var1_new = obj.BCR_s.aff.q1(2) + normrnd(0,MutStd);
                                %                                 q1_var2_new = obj.BCR_s.aff.q1(3) + normrnd(0,MutStd);
                                if(q1_var1_new<0) q1_var1_new=0; end
                                if(q1_var2_new<0) q1_var2_new=0; end
                                if( ( q1_var1_new >= 0) & ( q1_var2_new >= 0) )
                                    BCR_s_Kid.aff.q1 = obj.BCR_s.aff.q1;
                                    BCR_s_Kid.aff.q1(2) = q1_var1_new;
                                    BCR_s_Kid.aff.q1(3) = q1_var2_new;
                                    FlagMutOk = 1;
                                end
                            end
                        else
                            
                            q1_new = obj.BCR_s.aff.q1 + normrnd(0,MutStd);
                            if( q1_new >0)
                                BCR_s_Kid.aff.q1 = q1_new;
                                FlagMutOk = 1;
                            end
                        end
                        %                         BCR_s_Kid.aff.q1 = abs(obj.BCR_s.aff.q1 + normrnd(0,MutStd));
                    case 3
                        q2_new = obj.BCR_s.aff.q2 + normrnd(0,MutStd);
                        if( q2_new > 0)
                            BCR_s_Kid.aff.q2 = q2_new;
                            FlagMutOk = 1;
                        end
                    case 4
                        if(rand<obj.BCR_s.aff.qMutP)
                            new_rate = normrnd(0,MutStd);
                            q1_new = obj.BCR_s.aff.q1 + new_rate;
                            if( q1_new > 0)
                                BCR_s_Kid.aff.q1 = q1_new;
                                %                                 k_new =  0.25/q1_new;
                                %                                 E_new = log(1/k_new);
                                %                                 BCR_s_Kid.aff.EAbAg = E_new;
                                FlagMutOk = 1;
                            end
                        else
                            EAbAg_new = obj.BCR_s.aff.EAbAg + normrnd(0,MutStd);
                            if( EAbAg_new > 0)
                                BCR_s_Kid.aff.EAbAg = EAbAg_new;
                                FlagMutOk = 1;
                            end
                            %                             k = exp(-BCR_s_Kid.aff.EAbAg);
                            %                             k_new = k + new_rate;
                            %                             if( k_new > 0)
                            %                                 E_new = log(1/k_new);
                            %                                 BCR_s_Kid.aff.EAbAg = E_new;
                            % %                                 q1_new = 0.25/k_new;
                            % %                                 BCR_s_Kid.aff.q1 = q1_new;
                            %                                 FlagMutOk = 1;
                            %                             end
                        end
                    case 5
                        if(obj.BCR_s.aff.BindFlag==9)
                            if(rand<obj.BCR_s.aff.wcv(1)) %mutation in the conserved part
                                EAbAg_cons_new = obj.BCR_s.aff.EAbAg(1) + obj.LogNorm(MutStd);
                                if(EAbAg_cons_new<-1) EAbAg_cons_new=-1; end
                                %                                 if( EAbAg_cons_new >= 0)
                                BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                                %                                     BCR_s_Kid.aff.rcv = obj.ModifyRCV(obj.BCR_s.aff.EAbAg(1),EAbAg_cons_new,rcv_kid,1);
                                BCR_s_Kid.aff.EAbAg(1) = EAbAg_cons_new;
                                FlagMutOk = 1;
                                %                                 end
                            else %mutation in the variable part ; We take them to have a correlation of -2D such that increase for antigen
                                % a results in a decrease for b
                                %                                 dE = normrnd(0,MutStd);
                                %                                 EAbAg_var1_new = obj.BCR_s.aff.EAbAg(2) + dE;
                                %                                 EAbAg_var2_new = obj.BCR_s.aff.EAbAg(3) - dE;
                                EAbAg_var1_new = obj.BCR_s.aff.EAbAg(2) + obj.LogNorm(MutStd);
                                EAbAg_var2_new = obj.BCR_s.aff.EAbAg(3) + obj.LogNorm(MutStd);
                                if(EAbAg_var1_new<-1) EAbAg_var1_new=-1; end
                                if(EAbAg_var2_new<-1) EAbAg_var2_new=-1; end
                                %                                 if( (EAbAg_var1_new<0) & ((obj.BCR_s.aff.EAbAg(1)>0) | obj.BCR_s.aff.EAbAg(3)) )
                                %                                     EAbAg_var1_new=0;
                                %                                 end
                                %                                 if(EAbAg_var2_new<0) EAbAg_var2_new=0; end
                                %                                 if( ( EAbAg_var1_new >= 0) & ( EAbAg_var2_new >= 0) )
                                BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                                BCR_s_Kid.aff.EAbAg(2) = EAbAg_var1_new;
                                BCR_s_Kid.aff.EAbAg(3) = EAbAg_var2_new;
                                FlagMutOk = 1;
                                %                                 end
                            end
                        else
                            EAbAg_new = obj.BCR_s.aff.EAbAg + obj.LogNorm(MutStd);
                            if( EAbAg_new > 0)
                                BCR_s_Kid.aff.EAbAg = EAbAg_new;
                                FlagMutOk = 1;
                            end
                        end
                        
                        %%%%%%%%%
                        
                    case 6
                        
                        if(obj.BCR_s.aff.BindFlag==9)
                            if(rand<obj.BCR_s.aff.wcv(1)) %mutation in the conserved part
                                q1_cons_new = obj.BCR_s.aff.q1(1) + obj.LogNorm(MutStd);
                                if(q1_cons_new<0) q1_cons_new=0; end
                                if( q1_cons_new >= 0)
                                    BCR_s_Kid.aff.q1 = obj.BCR_s.aff.q1;
                                    BCR_s_Kid.aff.q1(1) = q1_cons_new;
                                    FlagMutOk = 1;
                                end
                            else %mutation in the variable part ; We take them to have a correlation of -2D such that increase for antigen
                                % a results in a decrease for b
                                %                                 dq = obj.LogNorm(MutStd);
                                %                                 q1_var1_new = obj.BCR_s.aff.q1(2) + dq;
                                %                                 q1_var2_new = obj.BCR_s.aff.q1(3) - dq;
                                q1_var1_new = obj.BCR_s.aff.q1(2) + obj.LogNorm(MutStd);
                                q1_var2_new = obj.BCR_s.aff.q1(3) + obj.LogNorm(MutStd);
                                if(q1_var1_new<0) q1_var1_new=0; end
                                if(q1_var2_new<0) q1_var2_new=0; end
                                if( ( q1_var1_new >= 0) & ( q1_var2_new >= 0) )
                                    BCR_s_Kid.aff.q1 = obj.BCR_s.aff.q1;
                                    BCR_s_Kid.aff.q1(2) = q1_var1_new;
                                    BCR_s_Kid.aff.q1(3) = q1_var2_new;
                                    FlagMutOk = 1;
                                end
                            end
                        else
                            
                            q1_new = obj.BCR_s.aff.q1 + obj.LogNorm(MutStd);
                            if( q1_new >0)
                                BCR_s_Kid.aff.q1 = q1_new;
                                FlagMutOk = 1;
                            end
                        end
                        
                        %%%%%%%%%
                        
                    case 7
                        if(obj.BCR_s.aff.BindFlag==9)
                            CorrVC = obj.BCR_s.aff.CorrVC; % the correlation between the effect of a mutation in the concerved site on the binding energy to the variable domain
                            CorrVV = obj.BCR_s.aff.CorrVV; % the correlation between the effect of a mutation in the variable site on the binding energy to the concerved domain
                            MutVC = obj.BCR_s.aff.MutVC; % the independent change in enrgy following a mutation in the concerved site on the binding energy to the variable domain
                            MutVV = obj.BCR_s.aff.MutVV; %

%                             Pvc = obj.BCR_s.aff.Pvc;
%                             Pvv = obj.BCR_s.aff.Pvv;
                            
                            dE = obj.LogNorm(MutStd);
                            Ec = obj.BCR_s.aff.EAbAg(1);
                            Ev1 = obj.BCR_s.aff.EAbAg(2);
                            Ev2 = obj.BCR_s.aff.EAbAg(3);
                            if(rand<obj.BCR_s.aff.wcv(1)) %mutation in the conserved part
                                
                                Ec = Ec + dE;
                                Ec = max(Ec,-1);
                                if(dE>0)
%                                     Ev1 = Ev1 -Pvc*dE + (1-Pvc)*obj.LogNorm(MutStd);
%                                     Ev2 = Ev2 -Pvc*dE + (1-Pvc)*obj.LogNorm(MutStd);
                                    Ev1 = Ev1 + CorrVC*dE + MutVC*obj.LogNorm(MutStd);
                                    Ev2 = Ev2 + CorrVC*dE + MutVC*obj.LogNorm(MutStd);
%                                     Ev1 = Ev1-dE;
%                                     Ev2 = Ev2-dE;
                                    Ev1 = max(Ev1,-1);
                                    Ev2 = max(Ev2,-1);
                                end
                                
                                BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                                BCR_s_Kid.aff.EAbAg(1) = Ec;
                                BCR_s_Kid.aff.EAbAg(2) = Ev1;
                                BCR_s_Kid.aff.EAbAg(3) = Ev2;
                                FlagMutOk = 1;
                                %                                 end
                            else %mutation in the variable part ; We take them to have a correlation of -2D such that increase for antigen
                                % a results in a decrease for b
                                %                                 dE = normrnd(0,MutStd);
                                %                                 EAbAg_var1_new = obj.BCR_s.aff.EAbAg(2) + dE;
                                %                                 EAbAg_var2_new = obj.BCR_s.aff.EAbAg(3) - dE;
   
                                % if the effect of affinity following
                                % mutation is negative, the mutation has
                                % equal change to reduce to the affinity
                                % torward each of the variable domains.
                                if(dE<0)
                                    if(rand<0.5)
                                        Ev1 = Ev1+dE;
                                    else
                                        Ev2 = Ev2+dE;
                                    end
                                else
                                    % if the effect of affinity following
                                % mutation is positive, the mutation has
                                % equal change to increase to the affinity
                                % torward each of the variable domains. If
                                % will reduce the affinity to the second
                                % variable domain and torward the conserved
                                % site.
%                                     Ec = Ec-dE;
%                                     Ec = Ec -Pvc*dE +(1-Pvc)*obj.LogNorm(MutStd);
                                    Ec = Ec + CorrVC*dE + MutVC*obj.LogNorm(MutStd);
                                                                       
                                    if(rand<0.5)
                                        Ev1 = Ev1+dE;
%                                         Ev2 = Ev2-dE;
%                                         Ev2 = Ev2-Pvv*dE + (1-Pvv)*obj.LogNorm(MutStd);
                                        Ev2 = Ev2 + MutVC*dE + MutVV*obj.LogNorm(MutStd);
                                    else
                                        Ev2 = Ev2+dE;
%                                         Ev1 = Ev1-dE;
%                                         Ev1 = Ev1-Pvv*dE + (1-Pvv)*obj.LogNorm(MutStd);
                                        Ev1 = Ev1 + MutVC*dE + MutVV*obj.LogNorm(MutStd);
                                    end
                                end
                                
                            end
                            Ec = max(Ec,-1);
                            Ev1 = max(Ev1,-1);
                            Ev2 = max(Ev2,-1);
                            BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                            BCR_s_Kid.aff.EAbAg(1) = Ec;
                            BCR_s_Kid.aff.EAbAg(2) = Ev1;
                            BCR_s_Kid.aff.EAbAg(3) = Ev2;
                            FlagMutOk = 1;
                            %                             end
                        else
                            EAbAg_new = obj.BCR_s.aff.EAbAg + obj.LogNorm(MutStd);
                            if( EAbAg_new > 0)
                                BCR_s_Kid.aff.EAbAg = EAbAg_new;
                                FlagMutOk = 1;
                            end
                        end
                    case 8
                        if(obj.BCR_s.aff.BindFlag==9)
                            CorrVC = obj.BCR_s.aff.CorrVC; % the correlation between the effect of a mutation in the concerved site on the binding energy to the variable domain
                            CorrVV = obj.BCR_s.aff.CorrVV; % the correlation between the effect of a mutation in the variable site on the binding energy to the concerved domain
                            MutVC = obj.BCR_s.aff.MutVC; % the independent change in enrgy following a mutation in the concerved site on the binding energy to the variable domain
                            MutVV = obj.BCR_s.aff.MutVV; %
%                             Pvc = obj.BCR_s.aff.Pvc;
%                             Pvv = obj.BCR_s.aff.Pvv;
                            
                            
                            Ec = obj.BCR_s.aff.EAbAg(1);
                            Ev1 = obj.BCR_s.aff.EAbAg(2);
                            Ev2 = obj.BCR_s.aff.EAbAg(3);
                            Ev3 = obj.BCR_s.aff.EAbAg(4);
                            
                            if(rand<obj.BCR_s.aff.wcv(1)) %mutation in the conserved part
%                                 dE = obj.LogNorm(MutStd);
                                dE = obj.LogNormCons(MutStd);
                                Ec = Ec + dE;
                                Ec = max(Ec,-1);
                                %                                 if(dE>0)
                                Ev1 = Ev1 + CorrVC*dE + MutVC*obj.LogNorm(MutStd);
                                Ev2 = Ev2 + CorrVC*dE + MutVC*obj.LogNorm(MutStd);
%                                 Ev1 = Ev1 -Pvc*dE + (1-Pvc)*obj.LogNorm(MutStd);
%                                 Ev2 = Ev2 -Pvc*dE + (1-Pvc)*obj.LogNorm(MutStd);
                                %                                     Ev1 = Ev1-dE;
                                %                                     Ev2 = Ev2-dE;
                                Ev1 = max(Ev1,-1);
                                Ev2 = max(Ev2,-1);
                                %                                 end
                                
                                BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                                BCR_s_Kid.aff.EAbAg(1) = Ec;
                                BCR_s_Kid.aff.EAbAg(2) = Ev1;
                                BCR_s_Kid.aff.EAbAg(3) = Ev2;
                                FlagMutOk = 1;
                                %                                 end
                            else %mutation in the variable part ; We take them to have a correlation of -2D such that increase for antigen
                                % a results in a decrease for b
                                %                                 dE = normrnd(0,MutStd);
                                %                                 EAbAg_var1_new = obj.BCR_s.aff.EAbAg(2) + dE;
                                %                                 EAbAg_var2_new = obj.BCR_s.aff.EAbAg(3) - dE;
   
                                % if the effect of affinity following
                                % mutation is negative, the mutation has
                                % equal change to reduce to the affinity
                                % torward each of the variable domains.
                                
%                                 if(dE<0)
%                                     if(rand<0.5)
%                                         Ev1 = Ev1+dE;
%                                     else
%                                         Ev2 = Ev2+dE;
%                                     end
%                                 else
                                    % if the effect of affinity following
                                % mutation is positive, the mutation has
                                % equal change to increase to the affinity
                                % torward each of the variable domains. If
                                % will reduce the affinity to the second
                                % variable domain and torward the conserved
                                % site.
%                                     Ec = Ec-dE;

             pvar=(copularnd('t',CorrVV,1,1));  %replace mutation scheme with copular to 

  %dE=zeros(length(pvar),1);
%dE(pvar<=0.0793848)=log(pvar(pvar<=0.0793848)/0.0793848)/1.43256;
%dE(pvar >0.0793848)=log(pvar((pvar>0.0793848))/0.920617)/-0.329412;
%dE=-dE;
             mu=1;                          %generate correlated uniform distributions
            sigma=.5;
            offset=4;
            dE = MutStd*(logninv(pvar,mu,sigma)-offset);  %take value of correlations then do inverse function to find value
            Ev1=Ev1+dE(1);
            Ev2=Ev2+dE(2);
            Ev3=Ev3+dE(3);

% % %                                     dE = obj.LogNorm(MutStd);
% % % 
% % %                                  
% % % %                                     Ec = Ec-Pvc*dE +(1-Pvc)*obj.LogNorm(MutStd);
% % % %                                     Ec = Ec-Pvc*dE +(1-Pvc)*obj.LogNormCons(MutStd);
% % % 
% % % global Virus;
% % %                                  if dE <0
% % % 
% % %                                     if min(min(obj.AffinityClass.aff.CorrVV)) ==1 && obj.AffinityClass.aff.MutVV ==0;
% % %                                        Ev1 = Ev1+dE;    %if on target or same keep them the same
% % %                                         Ev2 = Ev2+dE;
% % %                                         Ev3= Ev3+dE;
% % %                                     else
% % %                                      if Virus>0
% % %                                         Ev1 = Ev1+dE;
% % % 
% % %                                      else 
% % %                                         Ev2 = Ev2+dE;
% % % 
% % %                                      end
% % %                                     end
% % %                                  else
% % % %                                     Ec = Ec + CorrVC*dE + MutVC*obj.LogNorm(MutStd);
% % %                                     if Virus>0
% % %                                         Ev1 = Ev1+dE;
% % %                                         Ev2 = Ev2 + CorrVV(2,1)*dE + MutVV*obj.LogNorm(MutStd);
% % %                                         Ev3 = Ev3 + CorrVV(3,1)*dE + MutVV*obj.LogNorm(MutStd);
% % % 
% % %                                     else
% % %                                         Ev2 = Ev2+dE;
% % %                                         Ev1 = Ev1 + CorrVV(1,2)*dE + MutVV*obj.LogNorm(MutStd);
% % %                                         Ev3 = Ev3 + CorrVV(3,2)*dE + MutVV*obj.LogNorm(MutStd);
% % % 
% % %                                      
% % %                                     end
% % % 
% % %                                  end
% % % %                                 end
% % %                                 
                            end
                            Ec = max(Ec,-1);
                            Ev1 = max(Ev1,-2);
                            Ev2 = max(Ev2,-2);
                            Ev3 = max(Ev3, -2);
                            BCR_s_Kid.aff.EAbAg = obj.BCR_s.aff.EAbAg;
                            BCR_s_Kid.aff.EAbAg(1) = Ec;
                            BCR_s_Kid.aff.EAbAg(2) = Ev1;
                            BCR_s_Kid.aff.EAbAg(3) = Ev2;
                            BCR_s_Kid.aff.EAbAg(4) = Ev3;
                            FlagMutOk = 1;
                            %                             end
                        else
%                             EAbAg_new = obj.BCR_s.aff.EAbAg + obj.LogNorm(MutStd);
%                             % We do not need such a stringent condition
% %                             if( EAbAg_new > 0)
%                             if( EAbAg_new > -1)
%                                 BCR_s_Kid.aff.EAbAg = EAbAg_new;
%                                 FlagMutOk = 1;
%                             end
                            EAbAg_new = obj.BCR_s.aff.EAbAg;
                            EAbAg_new(2) = obj.BCR_s.aff.EAbAg(2) + obj.LogNorm(MutStd);
                            % We do not need such a stringent condition
%                             if( EAbAg_new > 0)
                            if( EAbAg_new > -1)
                                BCR_s_Kid.aff.EAbAg = EAbAg_new;
                                FlagMutOk = 1;
                            end
                        end
                        
                end
                if(FlagMutOk)
                    
                    Kid = Bcell(BCR_s_Kid);
                    Kid.BCR_s.ParentNum = obj.BCR_s.ParentNum+1;
                    obj.Sons = [obj.Sons ; Kid.SelfNum];
                else
                    Kid = [];
                end
            else
                Kid = [];
            end
        end
        
        %%%%%%%%%%%%%
        
        %%%%%%%%%%%%%
        
        function Kid = BirthJackpotKidMut(obj)
            
            X = obj.X;
            Y = int8(X'+1);
            cellJ = num2cell(Y);
            
            %             lambda = lambdaX(cellJ{:});
            lambda = obj.lambda;
            
            if(length(lambda)>1)
                Kid = [];
                return;
            end
            
            if(obj.DevLeft==0)
                obj.Flagln = 0;
            end
            
            if(obj.Flagln)
                lambda = obj.alpha*lambda;
            end
            if(  (lambda>obj.lambda0) && (obj.Flagln==0))
                obj.timeln=0;
                obj.DevLeft = obj.Dev0;
                obj.Flagln=1;
                lambda = obj.alpha*lambda;
            end
            
            p = obj.dt*lambda;
            if(rand <= p)
                MutStd = sqrt(2*obj.Diff*obj.dt);
                lambdaSon = abs(obj.lambda + normrnd(0,MutStd));
                %                 Kid = Bcell(lambdaSon,obj.mu,obj.Diff,obj.X,obj.dt,obj.SelfNum,obj.Jackpotflag,obj.lambda0,obj.Nmax,obj.lambdaBase,obj.BondsNum,obj.NBCR);
                BCR_s_Kid = obj.BCR_s;
                BCR_s_Kid.l = lambdaSon;
                Kid = Bcell(BCR_s_Kid);
                Kid.DevLeft = obj.DevLeft -1;
                if(Kid.DevLeft>0)
                    Kid.Flagln = 1;
                end
                obj.DevLeft = obj.DevLeft -1;
                obj.Sons = [obj.Sons ; Kid.SelfNum];
            else
                Kid = [];
            end
        end
        
        %%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%
        
        function Kid = BirthTXn(obj,lambdaX,ntot)
            
            Nmax = obj.Nmax;
            
            X = obj.X;
            
            Y = int8(X'+1);
            cellJ = num2cell(Y);
            lambda = lambdaX(cellJ{:});
            
            
            if(length(lambda)>1)
                Kid = [];
                return;
            end
            
            r = lambda - obj.mu;
            rn = r*(1- ntot/Nmax);
            if(rn<0)
                rn = 0;
            end
            p = obj.dt*(rn + obj.mu);
            
            if(rand <= p)
                Kid = BDD_moduleT(obj.lambda,obj.mu,obj.D,obj.X,obj.dt,obj.SelfNum);
                obj.Sons = [obj.Sons ; Kid.SelfNum];
            else
                Kid = [];
            end
        end
        
        
        %%%%%%%%%%%%%
        
        
        function Dflag = DeathT(obj);
            mu = obj.mu;
            x = obj.x;
            p = obj.dt*mu;
            if(rand <= p)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Dflag = DeathExpansionDecay(obj);
            mu = obj.BCR_s.muExpDecay;
            x = obj.x;
            p = obj.dt*mu;
            if(rand <= p)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Dflag = DeathTXn(obj,lambdaMeanAll,ntot);
            
            Nmax = obj.Nmax;
            
%             lambda = obj.lambdaBase;
            
            lambda = lambdaMeanAll;
            
            mu = obj.mu;
            
            r = lambda - mu;
            m_ntot = obj.mu + r*ntot/Nmax;
            p = obj.dt*m_ntot;
            if(rand <= p)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Dflag = DeathTInterCloneCom(obj,lambdaMeanAll,ntot,ni)
            
            Nmax = obj.Nmax;
            alpha = 0.1;
            
            ntotOther = ntot - ni;
            
            lambda = lambdaMeanAll;
            mu = obj.mu;
            
            r = lambda - mu;
            %             m_ntot = obj.mu + r*ntot/Nmax;
            m_ntot = obj.mu + r*(ni + alpha*ntotOther)/Nmax;
            p = obj.dt*m_ntot;
            if(rand <= p)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Dflag = DeathTXnTime(obj,lambdaX,ntot);
            
            Nmax = obj.Nmax;
            %             lambda = obj.lambda;
            X = obj.X;
            Y = int8(X'+1);
            cellJ = num2cell(Y);
            
            %             lambda = lambdaX(cellJ{:});
            lambda = obj.lambda;
            if(obj.Flagln)
                lambda = lambda*(1+exp(-obj.timeln/obj.tauNonlinear));
            end
            
            mu = obj.mu;
            
            r = lambda - mu;
            m_ntot = obj.mu + r*ntot/Nmax;
            p = obj.dt*m_ntot;
            if(rand <= p)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        
        
        function Dflag = DeathTn(obj,ntot);
            %             mu = obj.mu;
            dlm = obj.lambda - obj.mu;
            nrelax = 100;
            dmu = dlm*(1-exp(-ntot/nrelax));
            x = obj.x;
            p = obj.dt*(obj.mu + dmu);
            if(rand <= p)
                Dflag=1;
            else
                Dflag=0;
            end
        end
        function AgCapture(obj);
            
            %             AgCap = obj.AffinityClass.CAgCaptureFunc();
            %             obj.AgCapm2 = obj.AgCapm1;
            %             obj.AgCapm1 = obj.AgCap;
            %             obj.AgCap = AgCap;
            
            dt = 24*obj.dt;
            SurvivalTSteps = round(16/dt);
            
            obj.AgCapmArr(1:SurvivalTSteps-1) = obj.AgCapmArr(2:SurvivalTSteps);
            
            AgCap = obj.AffinityClass.CAgCaptureFunc();
            obj.AgCapmArr(SurvivalTSteps) = AgCap;
            obj.AgCap = AgCap;
            
        end
        
        function AgCaptureInit(obj);
            dt = 24*obj.dt;
            SurvivalTSteps = round(16/dt);
            for i=1:SurvivalTSteps
                obj.AgCapmArr(i) = obj.AffinityClass.CAgCaptureFunc();
            end
            obj.AgCap = obj.AgCapmArr(SurvivalTSteps);
        end
        
        function Dflag = AgCaptureThreshold(obj,AgThreshold);
            Dflag = 1;
            
            if(length(find(obj.AgCapmArr>=AgThreshold)))
                Dflag = 0;
            end
            
            %             tic
            %             for i=1:length(obj.AgCapmArr)
            %                 if(obj.AgCapmArr(i)>=AgThreshold)
            %                     Dflag = 0;
            %                 end
            %             end
            %             toc
        end
        
        % % poisson.m simulates a Poisson process
        %
        % lambda=1;      % arrival rate
        % Tmax=10;         % maximum time
        %
        % clear T;
        % T(1)=random('Exponential',1/lambda);
        % i=1;
        %
        % while T(i) < Tmax,
        %   T(i+1)=T(i)+random('Exponential',1/lambda);
        %   i=i+1;
        % end
        %
        % T(i)=Tmax;
        %
        % stairs(T(1:i), 0:(i-1));
        %
        %
        %             p = obj.lambda*obj.dt;
        %             if(rand
        %             A = obj.dt;
        %             noise_dev = sqrt(2*A);
        %             f = normrnd(0,noise_dev)'
        %             dx = f;
        %             obj.x = obj.x + dx;
        % %             coordnm = obj.coord.*repmat(obj.nm(1:size(obj.coord,2)),size(obj.coord,1),1);
        %         end
        %         function coordnm = px2nm(obj);
        %             coordnm = obj.coord.*repmat(obj.nm(1:size(obj.coord,2)),size(obj.coord,1),1);
        %         end
        %         function coordpx = nm2px(obj);
        %             coordpx = obj.coord./repmat(obj.nm,size(obj.coord,1),1);
        %         end
        function rcv = ModifyRCV(obj,Wold,Wnew,rcv,Agpart)
            dx=0.02;
            if(Wold<Wnew)
                if(Agpart==1)
                    rcv(1) = abs(rcv(1)+dx);
                    rcv(2) = abs(rcv(2)-dx);
                    rcv = rcv/sum(rcv);
                else
                    rcv(1) = abs(rcv(1)-dx);
                    rcv(2) = abs(rcv(2)+dx);
                    rcv = rcv/sum(rcv);
                end
            end
        end
        
        function a = LogNorm(obj,MutStd)
            mu     =  1; % lognormal mean
            sigma  =  0.5;  % lognormal standard deviation
            offsett = 4;
            a = MutStd*(lognrnd(mu,sigma,1)-offsett);
        end
        
        function a = LogNormCons(obj,MutStd)
            mu     =  1; % lognormal mean
            sigma  =  0.5;  % lognormal standard deviation
            offsett = 4;
            a = MutStd*(lognrnd(mu,sigma,1)-offsett);
        end
        
    end
end
