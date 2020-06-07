clear; clc;

ProjectPath = '/Volumes/Elham-Unifr/LongtermProject/ARCProject/DataAnonymization/EEGData/';

addpath(genpath('/Users/elhamb/Documents/Codes/Git/ARC'));
ResultsPath = '/Volumes/Elham-Unifr/LongtermProject/ARCProject/Results';
specanalysis = false;
pfanalysis = false;

%% Select the subset of subject for group level analysis
SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;34;37;38;44;47;48;52;61;65;66;70;72;78;]; %these are the subject with REC and HSMT and Plast: 49 and 63 -> good subjects without plast
%SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;37;38;44;47;48;61;65;66;70;72;78;];

SubIDs = arrayfun(@(x) ['ss-' num2str(x)],SubIDs,'uni',false);
load('Subjectinfo.mat');
SubIndex = find(ismember({SubjectData.SubID},SubIDs).*([SubjectData.Longitude]==0));

%% Spectrum estimation
if specanalysis
    EpochLen = 2500;
    Overlap = 1250;
    FreqBand = [1 15];
    ARC.Spectrumanalysis(ProjectPath ,'EpochLen' ,EpochLen ,'Overlap' ,Overlap , 'FreqBand', FreqBand,'RedoifExist',true,...
        'PlotResults',false,'Space','Source');
    
     ARC.Spectrumanalysis(ProjectPath ,'EpochLen' ,EpochLen ,'Overlap' ,Overlap , 'FreqBand', FreqBand,'RedoifExist',false,...
        'FigurePath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results','Space','Electrode');
end

%% Sensor space PARAFAC
% First Step: Calculate ARCs for REC conditions
if pfanalysis
    % SENSOR SPACE
    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','REC',...
        'Corcondia',false,'ResultsPath',ResultsPath,...
        'SubjectSelect',SubIndex);
    
    calculate ARCs for HSMT and Plast based on fixed frequency of REC condition loadings extracted from REC condition
    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','Plast',...
        'Corcondia',false,'ResultsPath',ResultsPath,'Space','Electrode',...
        'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');

    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','HSMT',...
        'Corcondia',false,'ResultsPath',ResultsPath,'Space','Electrode',...
    'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');
    
    % SOURCE SPACE
    % calculate ARCs for HSMT and Plast based on fixed
    % frequency of REC condition loadings extracted from REC condition in
    % source space
    
    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','REC',...
        'Corcondia',false,'ResultsPath',ResultsPath,'Space','Source',...
        'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');
    
    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','Plast',...
        'Corcondia',false,'ResultsPath',ResultsPath,'Space','Source',...
        'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');

    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','HSMT',...
        'Corcondia',false,'ResultsPath',ResultsPath,'Space','Source',...
    'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');
end

%% Group level analysis in sensor space
cond = 1;% 1= rec and hsmt, 2= rec plast hsmt
if cond ==1
    FileNames = {'REC1REC2REC3REC4','HSMTHSMT'};
    ModelNames = {'REC','HSMT'};
else
    FileNames = {'REC1REC2REC3REC4','Plast','HSMTHSMT'};% %
    ModelNames = {'REC','Plast','HSMT'};
end
%
Space = 'Electrode';
ModelPath = fullfile(ResultsPath,['PARAFAC_' Space]);
StatResults = ARC.ParafacAnova(ModelPath,SubjectData(SubIndex),'FileNames',FileNames,...
    'ModelNames',ModelNames,'ResultsPath',ResultsPath,'Space',Space,'PermNum',500,'redoAnalysis',false,'redoANOAVfigs',false,'plotasd',true);
%% Group level analysis in sensor space

Space = 'Source';
ModelPath = fullfile(ResultsPath,['PARAFAC_' Space]);
StatResults = ARC.ParafacAnova(ModelPath,SubjectData(SubIndex),'FileNames',FileNames,...
    'ModelNames',ModelNames,'ResultsPath',ResultsPath,'Space',Space,'PermNum',500,'redoAnalysis',false,'redoANOAVfigs',false,'plotasd',true);


