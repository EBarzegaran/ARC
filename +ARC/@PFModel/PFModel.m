classdef PFModel
    
    properties
        SubjectInfo
        Model
        Modes
        Freq
        Elabel
        Conditions
        CondLength
    end
    
    methods
        function obj = PFModel(SubjectInfo, Model,Modes,Freq,Elabel,Conditions,CondLength)
            % Initializing a PFModel object
            if exist('SubjectInfo','var')
                obj.SubjectInfo = SubjectInfo;
            else
                obj.SubjectInfo = [];
            end
            
            if exist('Model','var')
                obj.Model = Model;
            else
                obj.Model = [];
            end
            
            if exist('Freq','var')
                obj.Freq = Freq;
            else
                obj.Freq = [];
            end
            
            if exist('Elabel','var')
                obj.Elabel = Elabel;
            else
                obj.Elabel = [];
            end
            
            if exist('Conditions','var')
                obj.Conditions = Conditions;
            else
                obj.Conditions = [];
            end
            
            if exist('CondLength','var')
                obj.CondLength = CondLength;
            else
                obj.Conditions = [];
            end
            
            if exist('Modes','var')
                obj.Modes = Modes;
            else
                obj.Modes = [];
            end
        end
        
    end
end