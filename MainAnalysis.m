clear; clc;

ProjectPath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\DataAnonymization\EEGData\';
addpath(genpath('C:\Users\Elhamkhanom\Documents\Codes\Git\ARC\'));
ResultsPath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results';
specanalysis = false;
pfanalysis = false;

%% Select the subset of subject for group level analysis
%SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;34;37;38;44;47;48;52;61;65;66;70;72;78;]; %these are the subject with REC and HO and Plast
SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;37;38;44;47;48;61;65;66;70;72;78;]; %these are the subject with REC and HO and Plast
SubIDs = arrayfun(@(x) ['ss-' num2str(x)],SubIDs,'uni',false);
load('Subjectinfo.mat');
SubIndex = find(ismember({SubjectData.SubID},SubIDs).*([SubjectData.Longitude]==0));

%% Spectrum estimation
if specanalysis
    EpochLen = 2500;
    Overlap = 1250;
    FreqBand = [1 15];
    ARC.Spectrumanalysis(ProjectPath ,'EpochLen' ,EpochLen ,'Overlap' ,Overlap , 'FreqBand', FreqBand,'RedoifExist',true,...
         'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results');
end

%% Sensor space PARAFAC
% First Step: Calculate ARCs for REC conditions
if pfanalysis
    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','REC',...
        'Corcondia',false,'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results',...
        'SubjectSelect',SubIndex);

    % Second step: calculate ARCs for HSMT and Plast based on fixed frequency of REC condition loadings extracted from REC condition

    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','Plast',...
        'Corcondia',false,'ResultsPath',ResultsPath,...
        'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');


    ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','HSMT',...
        'Corcondia',false,'ResultsPath',ResultsPath,...
    'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');

end

%% Group level analysis

FileNames = {'REC1REC2REC3REC4','HSMTHSMT','Plast'};
ModelNames = {'REC','HSMT','Plast'};
ModelPath = fullfile(ResultsPath,['PARAFAC_' 'Electrode']);
ARC.ParafacAnova(ModelPath,SubjectData(SubIndex),'FileNames',FileNames,'ModelNames',ModelNames,'ResultsPath',ResultsPath);

