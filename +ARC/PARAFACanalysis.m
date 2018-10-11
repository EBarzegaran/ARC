function PARAFACanalysis(ProjectPath,varargin)

opt = ParseArgs(varargin,...
    'EpochLen'   ,2500,...
    'Overlap'    ,1250,...
    'FreqBand'  ,[5 15],...
    'Space'      ,'Electrode',...
    'Conditions' ,[],...
    'Subjectinfo'   ,[],...
    'SubjectSelect' ,[]...
        );
    
EpLen = opt.EpochLen; % Epoch length
MovWin = opt.Overlap;%Ovelpa of the moving window
addpath(genpath('C:\Users\ebarzega\Documents\My Works\TOOL\CSDtoolbox'));

if isempty(opt.Subjectinfo)
    load('+ARC\Private\Subjectinfo.mat');
    opt.Subjectinfo = SubjectData;
end

if isempty(opt.SubjectSelect)% select the first recordings of all subjects
    opt.SubjectSelect = find([SubjectData(:).Longitude]==0);
end

%% Read files and calculate spectrums and save them as RES class
if ~exist(fullfile(ProjectPath ,'FFTData'),'dir')
    mkdir(fullfile(ProjectPath ,'FFTData'));
end
for S = 1:numel(opt.SubjectSelect)
    sub = opt.SubjectSelect(S);
    display([SubjectData(sub).SubID]);
    temp = ARC.RES();
    if ~temp.loadRES(fullfile(ProjectPath ,'FFTData'),SubjectData(sub))
        [RESdata] = RawEEGtoRES(ProjectPath,SubjectData(sub),opt.EpochLen ,opt.Overlap ,opt.FreqBand(2),opt.FreqBand(1));
        RESdata.saveRES(fullfile(ProjectPath ,'FFTData'));
    end
end

end