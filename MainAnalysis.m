clear; clc;

ProjectPath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\DataAnonymization\EEGData\';
addpath(genpath('C:\Users\Elhamkhanom\Documents\Codes\Git\AlphaComponentAnalysis\'));
%% Analysis Parameters
% Epoching parameters
EpochLen = 2500;
Overlap = 1250;
FreqBand = [1 15];

%% Spectrum estimation
ARC.Spectrumanalysis(ProjectPath ,'EpochLen' ,EpochLen ,'Overlap' ,Overlap , 'FreqBand', FreqBand,...
    'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results');

%% Sensor space PARAFAC
ARC.PARAFACanalysis(ProjectPath ,'EpochLen' ,EpochLen ,'Overlap' ,Overlap , 'FreqBand', [5 15],'VarianceMode','temporal',...
    'Corcondia',false,'ResultsPath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results');

