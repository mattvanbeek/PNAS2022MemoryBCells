classdef ImmunogenicityAgClass < handle
    properties
        AgThGCSeedHead; % required amount of Ag need to enter a GC when the Ag contains the HA head
        AgThGCSeedStemOnly; % required amount of Ag need to enter a GC when the Ag does not contain the HA head
    end
    methods
        function obj = ImmunogenicityAgClass(AgThGCSeedHead,AgThGCSeedStemOnly)
            obj.AgThGCSeedHead = AgThGCSeedHead;
            obj.AgThGCSeedStemOnly = AgThGCSeedStemOnly;
        end
        function AgThresholdGCSeed = ImmunogenAgTh(obj,Ag)
            switch Ag
                case '1HA'
                    AgThresholdGCSeed = obj.AgThGCSeedHead;
                case '40VirusHA'
                    AgThresholdGCSeed = obj.AgThGCSeedHead;
                case 'NPHA'
                    AgThresholdGCSeed = obj.AgThGCSeedHead;
                case 'NPSS'
                    AgThresholdGCSeed = obj.AgThGCSeedStemOnly;
            end
            AgThresholdGCSeed=100;
        end
    end
end