classdef PFModel
    % a class to store the output of PARAFAC analysis
    properties
        SubjectInfo
        Model
        Modes
        VarMode
        Freq
        Elabel
        Conditions
        CondLength
    end
    properties (Dependent)
       ARpeaks
       ARwidths
    end
    
    methods
        %% creat, save and load functions
        function obj = PFModel(SubjectInfo, Model,Modes,VarMode,Freq,Elabel,Conditions,CondLength)
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
            if exist('VarMode','var')
                obj.VarMode = VarMode;
            else
                obj.VarMode = [];
            end
        end
        
        function Name = SaveName(obj)
            Name = ['PARAFAC_' obj.SubjectInfo.SubID '_RecordNum' num2str(obj.SubjectInfo.Longitude+1) '_' [obj.Conditions{:}] '_' obj.VarMode];
        end
        
        function savePFModel(obj,Path)
            % Saves the object in the deifined path
            try
                SaveName = obj.SaveName;
            catch
                error('SubjectInfo property should have SubID and Longitude fields');
            end
            try
                if ~exist('Path','var')
                    Path = [];
                end
                save(fullfile(Path,SaveName),'obj');
            catch
            end
        end
        
        function obj = loadPFModel(obj,Path, SubjectInfo,Conditions,VarMode)
            %loads the RES of the subject from the defined folder
            if ~exist('Conditions','var')
                Conditions = '*';
            else
                if iscell(Conditions)
                    Conditions = [Conditions{:}];
                end
                if exist('VarMode','var')
                    Conditions = [Conditions '_' VarMode];
                end
            end
            try
                SaveName = ['PARAFAC_' SubjectInfo.SubID '_RecordNum' num2str(SubjectInfo.Longitude+1) '_'];
            catch
                error('SubjectInfo should have SubID and Longitude fields');
            end
                if ~exist('Path','var')
                    Path = [];
                end
                File = subfiles([fullfile(Path,[SaveName Conditions]) '*.mat'],1);
                if ischar(File{1})
                    if numel(File)>1
                        warning('More than one PARAFAC model exist for this subject. Please indicate the conditions and variance mode')
                    end
                    load(File{1});
                else
                    obj=[];
                end
        end

        %% Dependent variables
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
                [m] = mean(Fmode); % OR min
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
        
        %% other functions
        function [Mean,SEM] = PFTempMeans(obj)
            TempLen = [0 cumsum(obj.CondLength)];
            Rang = [TempLen(1:end-1)+1;TempLen(2:end)];
            Temp= obj.Model{(strcmpi(obj.Modes,'temporal'))};
            TempM = arrayfun(@(x) mean(Temp(Rang(1,x):Rang(2,x),:)),1:size(Rang,2),'uni',false);
            TempS = arrayfun(@(x) std(Temp(Rang(1,x):Rang(2,x),:)),1:size(Rang,2),'uni',false);
            Mean = cat(1,TempM{:});
            SEM = cat(1,TempS{:})./repmat(sqrt(diff(Rang)+1)',[1 size(TempS{1},2)]);
        end
        
        function obj = OrganizeARCs(obj)
            % order according to peak
            % if the peak is <6Hz or >14Hz, remove that component
            % if the width is zeros, remove that component
            % THIS FUNCTION SHOULD BE UPDATED LATER
            NumCom = obj.SubjectInfo.Compnum;
            Comp = obj.SubjectInfo.Comp;
            Comp = Comp(Comp>0);
            Peaks = obj.ARpeaks(Comp);
         
            % order the components according to their frequency peaks
            [~,Order] = sort(Peaks,'descend');
            if numel(Order)>2
                Order(2:3) = Order(3:-1:2);
            end
            Comp = Comp(Order);
            obj.SubjectInfo.Comp = Comp;
        end
        
        function Loading = GetLoading(obj,Comp,Mode)
            % Returns the loading of component #Comp with the Mode
            % Mode: 'spatial'/['frequency']/'temporal';
            if ~exist('Mode','var')
                Mode = 'frequency';
            end
            obj = obj.OrganizeARCs;
            objComp = obj.SubjectInfo.Comp; objComp(objComp==0)=[];
            Loading = obj.Model{strcmpi(obj.Modes,Mode)}(:,objComp);
            Loading = Loading(:,Comp);
        end   
    end
end