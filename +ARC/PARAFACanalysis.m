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
    'VarianceMode'  ,'spatial'...
        );

if isempty(opt.Subjectinfo)
    load('+ARC\Private\Subjectinfo.mat');
    opt.Subjectinfo = SubjectData;
end

if isempty(opt.SubjectSelect)% select the first recordings of all subjects
    opt.SubjectSelect = find([SubjectData(:).Longitude]==0);
end

%% Read RES class
for S = 1:numel(opt.SubjectSelect)
    sub = opt.SubjectSelect(S);
    display([SubjectData(sub).SubID]);
    temp = ARC.RES();
    if isempty(temp.loadRES(fullfile(ProjectPath ,'FFTData'),SubjectData(sub)))
        error('RES class for subject not found, please run ARC.Spectrumanalysis first');
    else
        RESdata{S} = temp.loadRES(fullfile(ProjectPath ,'FFTData'),SubjectData(sub));
    end

    Model = ResParafac(RESdata{S},'FreqBand',opt.FreqBand,'Conditions',opt.Conditions,...
        'Electrodes',opt.Electrodes,'Corconia',opt.Corcondia,'VarianceMode',opt.VarianceMode);
end
end