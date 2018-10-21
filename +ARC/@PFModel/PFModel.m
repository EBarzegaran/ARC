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
    properties (Dependent)
       ARpeaks
       ARwidths
       ARTempMeans
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
        
        function value = get.ARpeaks(obj)
            % find the Modes which is frequency
            M = find(strcmpi(obj.Modes,'frequency'));
            % Then find their peak... if it is minimum or maximum?
            if ~isempty(M)
                Fmode = obj.Model{M};
                [~,Ind] = max(Fmode);
                value = obj.Freq(Ind);
            else
                value = [];
            end
        end
        
        function value = get.ARwidths(obj)
            % find the Modes which is frequency
            M = find(strcmpi(obj.Modes,'frequency'));
            % Then find full width half maximym...
            if ~isempty(M)
                Fmode = obj.Model{M};
                [M,MI] = max(Fmode);
                [m] = mean(Fmode); % OR mean
                HM= (M+m)/2; % half maximum
                for c = 1:numel(HM)
                    ind1 = diff(Fmode(:,c)<HM(c)); 
                    I1=find(ind1==-1)'; I2=find(ind1==1)';
                    if numel(I1)==numel(I2)
                        I = [I1;I2];
                        % fin the peak with global maximum
                        Ind = MI(c)>=I(1,:).*MI(c)<=I(2,:);
                        Width(c) = (I2(Ind)-I1(Ind)).*diff(obj.Freq(1:2));
                    else
                        Peak(c) = 1;
                        Width(c) = 0;
                    end 
                end
                value = Width;
            else
                value = [];
            end
        end
        
        function value = get.ARTempMeans(obj)
            TempLen = [0 cumsum(obj.CondLength)];
            Rang = [TempLen(1:end-1)+1;TempLen(2:end)];
            Temp= obj.Model{(strcmpi(obj.Modes,'temporal'))};
            TempM = arrayfun(@(x) mean(Temp(Rang(1,x):Rang(2,x),:)),1:size(Rang,2),'uni',false);
            value = cat(1,TempM{:});
        end
    end
end