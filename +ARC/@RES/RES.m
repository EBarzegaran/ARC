classdef RES
    % Resting-state EEG Spectrum
    properties
        SubjectInfo
        FFTData
        Freq
        EpochLen
        Elabels
        Overlap
    end
    
    methods
        %% creat, save and load functions
        function obj = RES(SubjectInfo, FFTData,CondNames, EpochLen, Overlap,Freq,Elabels)
            % initialize a RES class object
            if exist('SubjectInfo','var')
                obj.SubjectInfo = SubjectInfo;
            else
                obj.SubjectInfo = [];
            end
            
            if exist('FFTData','var')
                obj.FFTData = cellfun(@(x,y) ARC.FData(x,y),FFTData,CondNames); % ARRAY of FData blocks
            else
                obj.FFTData = [];
            end
            
            if exist('Freq','var')
                obj.Freq = Freq;
            else
                obj.Freq = [];
            end
            
            if exist('EpochLen','var')
                obj.EpochLen = EpochLen;
            else
                obj.EpochLen = [];
            end
            
            if exist('Elabels','var')
                obj.Elabels = Elabels;
            else
                obj.Elabels = [];
            end
            
            if exist('Overlap','var')
                obj.Overlap = Overlap;
            else
                obj.Overlap = [];
            end 
        end
        
        function saveRES(obj,Path)
            % Saves the object in the deifined path
            try
                SaveName = ['FFTData_' obj.SubjectInfo.SubID '_RecordNum' num2str(obj.SubjectInfo.Longitude+1)];
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
        
        function obj = loadRES(obj,Path, subinfo)
            %loads the RES of the subject from the defined folder
            try
                SaveName = ['FFTData_' subinfo.SubID '_RecordNum' num2str(subinfo.Longitude+1) '.mat'];
            catch
                error('SubjectInfo should have SubID and Longitude fields');
            end
                if ~exist('Path','var')
                    Path = [];
                end
                if exist(fullfile(Path,SaveName),'file')
                    load(fullfile(Path,SaveName));
                else
                    obj=[];
                end

        end
        
        %% Visualization function
        function FigHandler = PlotSpectrum(obj,varargin)
            % INPUTS (optional):
            % Varargin:
            % FreqBand: is a vector [LF HF], where LF is the minimum and HF
                         % is maximum frequency of the plots
            % Electrodes: Index of Electrodes to be plotted, each of them
                        % is plotted in a 2 x numeber of electrodes/2 plot. Better for
                        % visulization of left and right electrodes
            % Conditions: can be any of the following options or a cell array containing a list of them:
                         % 'REC': average of all resting eyes close conditions
                         % 'REC1':  first resting eyes closed
                         % 'REC2' or 'REC-Plast': REC after plastisine
                         % 'REC3' or 'REC-HSMT': REC after haptic comparison
                         % 'REC4' or 'REC-HN': REC after haptic navigation
                         % 'REO': rest eyes open
                         % 'Plast': Plastisine condition
                         % 'HSMT': Haptic object matching task
                         % 'HN': Haptic Navigation
            % Colors: a numel(Condition) x 3 array indicating colors of
                        % each condition
            % Normalize: if 1, normalize yaxis of all electrode spectrums to maximum
            % value in these electrodes
            opt = ParseArgs(varargin,...
                'Normalize'     ,1,...
                'FreqBand'      ,[min(obj.Freq) max(obj.Freq)],...
                'Electrodes'    ,[15 17 24 28 25 27 54 57 30 31],...
                'Conditions'    ,{'REC','Plast','HSMT'},...
                'Colors'        ,[],...
                'SavePath'      ,[]...
                );
            
            % Condition list and default conditions
            CondList = {'REC','REC1','REC2','REC-Plast','REC3','REC-HSMT','REC4','REC-HN','REO','Plast','HSMT','HN'};
            if ~prod(ismember(opt.Conditions,CondList))
                error('Undefined Conditions for plotting');
            end

            
            % Set default colors
            defaultcolors = struct('Conditions',CondList,...
                    'Colors',num2cell([0 0 0; .1 .1 .1; .3 .3 .3; .3 .3 .3; .5 .5 .5 ; .5 .5 .5; .7 .7 .7 ; .7 .7 .7; 0 1 0 ;1 0 1; 0 .3 1; 1 0 0],2)');
            if isempty(opt.Colors)
                  I = cellfun(@(x) find(strcmp(CondList,x)),opt.Conditions,'uni',false);           
                  opt.Colors = {defaultcolors(cell2mat(I)).Colors};
                  opt.Colors = cat(1,opt.Colors{:});
            end
            
            % prepare the power spectrum for plotting
            objConds = {obj.FFTData(:).Condition};
            for c = 1:numel(opt.Conditions)
                 switch opt.Conditions{c}
                     case 'REC' 
                         Inds{c} = find(contains(lower(objConds),'rec'));
                     case 'REC1'
                         Inds{c} = find(strcmpi(objConds,'rec'));
                     case {'REC2','REC-Plast'}
                         Inds{c} = find(strcmpi(objConds,'REC-mu'));
                     case {'REC3','REC-HSMT'}
                         Inds{c} = find(strcmpi(objConds,'REC-hob'));
                     case {'REC4','REC-HN'}
                         Inds{c} = find(strcmpi(objConds,'REC-hn'));
                     case 'REO'
                         Inds{c} = find(strcmpi(objConds,'REO'));
                     case 'Plast'
                         Inds{c} = find(strcmpi(objConds,'Plast'));
                     case 'HSMT'
                         Inds{c} = find(contains(lower(objConds),{'hob1','hob2'}));
                     case 'HN'
                         Inds{c} = find(contains(lower(objConds),{'hn1','hn2'}));
                 end
                 
                 S = {obj.FFTData(Inds{c}).Data};
                 Spec{c} = mean(sqrt(abs(cat(3,S{:}))),3);
            end
            
            if opt.Normalize
                Specall = cat(3,Spec{:});
                Ind =(obj.Freq>=opt.FreqBand(1)).*(obj.Freq<=opt.FreqBand(2));
                Specall = Specall(:,Ind==1,:);
                MYaxis = max(Specall(:))*1.1;
            end
            
            %Plot the spectrum
            FS = 14;
            FigHandler = figure;
            ElNum = floor(numel(opt.Electrodes)/2)*2;
            for sp = 1:ElNum
                subplot(ElNum/2,2,sp);
                title(obj.Elabels{opt.Electrodes(sp)},'fontsize',FS);
                % NOt empty spec
                hold on;
                for c = 1:numel(opt.Conditions)
                    if ~isempty(Spec{c})
                        plot(obj.Freq,squeeze(Spec{c}(opt.Electrodes(sp),:)),'linewidth',2,'color',opt.Colors(c,:));
                    end
                end
                xlim(opt.FreqBand);
                if opt.Normalize
                    ylim([0 MYaxis]);
                end
                if sp == 2
                    L = legend(opt.Conditions);
                    set(L,'Position',get(L,'Position')+[.05 .05 0 0]);
                end
                
                if sp==(floor(ElNum/4)*2+1)
                    ylabel(['ASD (\muv / ' 'Hz)'],'fontsize',FS);
                end
                
                if sp ==ElNum
                    xlabel('Frequency(Hz)','fontsize',FS);
                end
            end
            
            if ~isempty(opt.SavePath)
                set(FigHandler,'PaperPosition',[1 1 7 ElNum]);
                EName = [obj.Elabels{opt.Electrodes}];
                CName = [opt.Conditions{:}];
                print(fullfile(opt.SavePath,['EEGSPectrum_' obj.SubjectInfo.SubID '_RecordNum' num2str(obj.SubjectInfo.Longitude+1) '_' EName '_' CName '.tif']),'-r300','-dtiff');
            end
        end
    end
end

