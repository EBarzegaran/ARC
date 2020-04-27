function PARAFACanalysis(ProjectPath,varargin)
% Apply PARAFAC on the RES data ..

opt = ParseArgs(varargin,...
    'FreqBand'      ,[5 15],...
    'Space'         ,'Electrode',...
    'Conditions'    ,[],...
    'Electrodes'    ,[],...
    'Corcondia'     ,false,...
    'Subjectinfo'   ,[],...
    'SubjectSelect' ,[],...
    'ResultsPath'   ,[],...
    'VarianceMode'  ,'spatial',...
    'FixedFreqLoading', false,...
    'FixedModel'      ,[],...
    'SaveFigures'   ,false...
        );

if isempty(opt.Subjectinfo)
    load('+ARC\Private\Subjectinfo.mat');
    opt.Subjectinfo = SubjectData;
end

if isempty(opt.SubjectSelect)% select the first recordings of all subjects
    opt.SubjectSelect = find([SubjectData(:).Longitude]==0);
end

%% Read RES class
if ~exist(fullfile(opt.ResultsPath,['PARAFAC_' opt.Space]),'dir')
    mkdir(fullfile(opt.ResultsPath,['PARAFAC_' opt.Space]));
end
if ~exist(fullfile(opt.ResultsPath,['PARAFAC_' opt.Space],'Figures'),'dir')
    mkdir(fullfile(opt.ResultsPath,['PARAFAC_' opt.Space],'Figures'));
end
for S = 1:numel(opt.SubjectSelect)
    % Read the subject's amplitude spectrum density, if already computed,
    % just loads them
    sub = opt.SubjectSelect(S);
    display([SubjectData(sub).SubID]);
    temp = ARC.RES();
    if isempty(temp.loadRES(fullfile(ProjectPath ,'FFTData'),SubjectData(sub),opt.Space))
        error('RES class for subject not found, please run ARC.Spectrumanalysis first');
    else
        RESdata{S} = temp.loadRES(fullfile(ProjectPath ,'FFTData'),SubjectData(sub),opt.Space);
    end

    % if PARAFAC should be done with fixed loading
    if opt.FixedFreqLoading
        ModelTemp = ARC.PFModel();
        TConds = RESdata{1}.FindConditions(opt.FixedModel);
        ModelTemp = ModelTemp.loadPFModel(fullfile(opt.ResultsPath,['PARAFAC_' 'Electrode']), RESdata{S}.SubjectInfo,RESdata{1}.GetCondNames([TConds{:}]),opt.VarianceMode,'Electrode');
    else
        ModelTemp=[];
    end
    Model = ResParafac(RESdata{S},'FreqBand',opt.FreqBand,'Conditions',opt.Conditions,...
        'Electrodes',opt.Electrodes,'Corconia',opt.Corcondia,'VarianceMode',opt.VarianceMode,...
        'FixedFreqLoading',opt.FixedFreqLoading,'FixedModel',ModelTemp,'Space',opt.Space);
    Model.savePFModel(fullfile(opt.ResultsPath,['PARAFAC_' opt.Space]));
    if opt.SaveFigures
        Model.plotPFmodel(1,fullfile(opt.ResultsPath,['PARAFAC_' opt.Space],'Figures'));
        close;
    end
end
end