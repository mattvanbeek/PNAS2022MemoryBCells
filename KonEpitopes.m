classdef KonEpitopes < handle
    properties
        NArr0;
        P1Arm0;
        P2Arm0;
        kOn1Arm0;
        kOn2Arm0;
        
        NMax;
        
        NArr;
        P1Arm;
        P2Arm;
        kOn1Arm;
        kOn2Arm;
        
        EpitopeNum;
        q01Arm;
        q02Arm;
    end
    methods
        function obj = KonEpitopes(NP,q01Arm,q02Arm);
            
            obj.q01Arm = q01Arm;
            obj.q02Arm = q02Arm;
            
            %             if(nargin > 0)
            if(strcmp(NP,'NPHA'))
                load('Flu_NPHA_082418');
%                 q0 = 0.2;
                
                obj.NArr0 = [100];
                
                obj.NMax = obj.NArr0(end);
                
                % The probability of one arm to be bound at the head
                % and the stem epitopes
                
                %                 % The probabilities of one arm to be bound at the head
                %                 % and the stem epitopes are equal
                obj.P1Arm0 = P1Arm_HA';
                obj.P2Arm0 = P2Arm_HA';
                
                obj.EpitopeNum = length(P1Arm_HA);
                
                % The on rate of one arm to be bound at the head
                % and the stem epitopes
                
                
                % The on rates of one arm to be arrive at the head
                % and the stem epitopes are equal
                
                obj.kOn1Arm0 = kOn1Arm_HA';
                obj.kOn2Arm0 = kOn2Arm_HA';
                
                % Giving the stem clone a fitness adventage to capture Ag
                % with the second Arm.
                
                obj.kOn1Arm0 = obj.q01Arm*obj.kOn1Arm0;
                obj.kOn2Arm0 = obj.q02Arm*obj.kOn2Arm0;
                
                obj.NArr0 = [ 0 obj.NArr0];
                
                obj.P1Arm0 = [  zeros(obj.EpitopeNum,1) obj.P1Arm0];
                obj.P2Arm0 = [ zeros(obj.EpitopeNum,1) obj.P2Arm0];
                obj.kOn1Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn1Arm0];
                obj.kOn2Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn2Arm0];
                obj.CalcOnRates();
            elseif(strcmp(NP,'NPSS'))
                load('Flu_NPSS_082118');
                
%                 q0 = 0.2;
                
                obj.NArr0 = [100];
                
                obj.NMax = obj.NArr0(end);
                
                % The probability of one arm to be bound at the head
                % and the stem epitopes
                
                obj.P1Arm0 = P1Arm_SS';
                obj.P2Arm0 = P2Arm_SS';
                
                obj.EpitopeNum = length(P1Arm_SS);
                
                
                % The on rate of one arm to be bound at the head
                % and the stem epitopes
                
                obj.kOn1Arm0 = kOn1Arm_SS';
                obj.kOn2Arm0 = kOn2Arm_SS';
                
                % Giving the stem clone a fitness adventage to capture Ag
                % with the second Arm.
                
                obj.kOn1Arm0 = obj.q01Arm*obj.kOn1Arm0;
                obj.kOn2Arm0 = obj.q02Arm*obj.kOn2Arm0;
                
                obj.NArr0 = [ 0 obj.NArr0];
                
                obj.P1Arm0 = [  zeros(obj.EpitopeNum,1) obj.P1Arm0];
                obj.P2Arm0 = [ zeros(obj.EpitopeNum,1) obj.P2Arm0];
                obj.kOn1Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn1Arm0];
                obj.kOn2Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn2Arm0];
                obj.CalcOnRates();
            elseif(strcmp(NP,'1HA'))
                load('Flu_1HA_041119');
%                 q0 = 0.2;
                
                obj.NArr0 = [100];
                
                obj.NMax = obj.NArr0(end);
                
                % The probability of one arm to be bound at the head
                % and the stem epitopes
                
                %                 % The probabilities of one arm to be bound at the head
                %                 % and the stem epitopes are equal
                obj.P1Arm0 = P1Arm_HA';
                obj.P2Arm0 = P2Arm_HA';
                
                obj.EpitopeNum = length(P1Arm_HA);
                
                % The on rate of one arm to be bound at the head
                % and the stem epitopes
                
                
                % The on rates of one arm to be arrive at the head
                % and the stem epitopes are equal
                
                obj.kOn1Arm0 = kOn1Arm_HA';
                obj.kOn2Arm0 = kOn2Arm_HA';
                
                % Giving the stem clone a fitness adventage to capture Ag
                % with the second Arm.
                
                obj.kOn1Arm0 = obj.q01Arm*obj.kOn1Arm0;
                obj.kOn2Arm0 = obj.q02Arm*obj.kOn2Arm0;
                
                obj.NArr0 = [ 0 obj.NArr0];
                
                obj.P1Arm0 = [  zeros(obj.EpitopeNum,1) obj.P1Arm0];
                obj.P2Arm0 = [ zeros(obj.EpitopeNum,1) obj.P2Arm0];
                obj.kOn1Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn1Arm0];
                obj.kOn2Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn2Arm0];
                obj.CalcOnRates();
                
                
                elseif(strcmp(NP,'40VirusHA'))
                %load('virusha40Spikes_081419_OffStemModified');
                %load('virusha50Spikes_081719');
                %load('virushaHead50Stem56_081819');
                load('virusha56Spikes_082319');
%                 q0 = 0.2;
                
                obj.NArr0 = [100];
                
                obj.NMax = obj.NArr0(end);
                
                % The probability of one arm to be bound at the head
                % and the stem epitopes
                
                %                 % The probabilities of one arm to be bound at the head
                %                 % and the stem epitopes are equal
                obj.P1Arm0 = P1Arm_HA';
                obj.P2Arm0 = P2Arm_HA';
                
                obj.EpitopeNum = length(P1Arm_HA);
                
                % The on rate of one arm to be bound at the head
                % and the stem epitopes
                
                
                % The on rates of one arm to be arrive at the head
                % and the stem epitopes are equal
                
                obj.kOn1Arm0 = kOn1Arm_HA';
                obj.kOn2Arm0 = kOn2Arm_HA';
                
                % Giving the stem clone a fitness adventage to capture Ag
                % with the second Arm.
                
                obj.kOn1Arm0 = obj.q01Arm*obj.kOn1Arm0;
                obj.kOn2Arm0 = obj.q02Arm*obj.kOn2Arm0;
                
                obj.NArr0 = [ 0 obj.NArr0];
                
                obj.P1Arm0 = [  zeros(obj.EpitopeNum,1) obj.P1Arm0];
                obj.P2Arm0 = [ zeros(obj.EpitopeNum,1) obj.P2Arm0];
                obj.kOn1Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn1Arm0];
                obj.kOn2Arm0 = [ zeros(obj.EpitopeNum,1) obj.kOn2Arm0];
                obj.CalcOnRates();
                
            end
            
        end
        
        function CalcOnRates(obj);
            
            NArr = [0:1:obj.NMax];
            obj.P1Arm = zeros(obj.EpitopeNum,obj.NMax+1);
            obj.P2Arm = zeros(obj.EpitopeNum,obj.NMax+1);
            obj.kOn1Arm = zeros(obj.EpitopeNum,obj.NMax+1);
            obj.kOn2Arm = zeros(obj.EpitopeNum,obj.NMax+1);
            for N=NArr
                for j=1:obj.EpitopeNum
                    obj.P1Arm(j,N+1) = interp1(obj.NArr0,obj.P1Arm0(j,:),N);
                    obj.P2Arm(j,N+1) = interp1(obj.NArr0,obj.P2Arm0(j,:),N);
                    obj.kOn1Arm(j,N+1) = interp1(obj.NArr0,obj.kOn1Arm0(j,:),N);
                    obj.kOn2Arm(j,N+1) = interp1(obj.NArr0,obj.kOn2Arm0(j,:),N);
                end
            end
        end
        
        %%%
        %%%
    end
    
end

