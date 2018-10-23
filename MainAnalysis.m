clear; clc;

ProjectPath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\DataAnonymization\EEGData\';
addpath(genpath('C:\Users\Elhamkhanom\Documents\Codes\Git\AlphaComponentAnalysis\'));
%% Analysis Parameters
% Epoching parameters
EpochLen = 2500;
Overlap = 1250;
FreqBand = [1 15];

%% Spectrum estimation
ARC.Spectrumanalysis(ProjectPath ,'EpochLen' ,EpochLen ,'Overlap' ,Overlap , 'FreqBand', FreqBand,'RedoifExist',true,...
     'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results');

%% Select the subset of subject for group level analysis
%SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;34;37;38;44;47;48;52;61;65;66;70;72;78;]; %these are the subject with REC and HO and Plast
SubIDs= [1;3;4;8;9;12;13;16;17;19;21;22;23;27;29;32;37;38;44;47;48;61;65;66;70;72;78;]; %these are the subject with REC and HO and Plast
SubIDs = arrayfun(@(x) ['ss-' num2str(x)],SubIDs,'uni',false);
load('Subjectinfo.mat');
SubIndex = find(ismember({SubjectData.SubID},SubIDs).*([SubjectData.Longitude]==0));

%% Sensor space PARAFAC
% First Step: Calculate ARCs for REC conditions

ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','REC',...
    'Corcondia',false,'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results',...
    'SubjectSelect',SubIndex);

%% Second step: calculate ARCs for HSMT and Plast based on fixed frequency
%loadings extracted from REC condition

ARC.PARAFACanalysis(ProjectPath ,'FreqBand', [5 15],'VarianceMode','temporal','SaveFigures',true,'Conditions','Plast',...
    'Corcondia',false,'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results',...
    'SubjectSelect',SubIndex, 'FixedFreqLoading', true, 'FixedModel','REC');
