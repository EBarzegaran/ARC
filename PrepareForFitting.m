clear; clc;

ProjectPath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\DataAnonymization\EEGData\';
addpath(genpath('C:\Users\Elhamkhanom\Documents\Codes\Git\ARC\'));
ResultsPath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results';
SavePath = 'C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results\FFTData_NNRfitting';
specanalysis = false;
pfanalysis = false;

%% Select subjects
load('Subjectinfo.mat');
SubIndex = find([SubjectData(:).Longitude]==0);
%% Spectrum estimation
SpecAnalysis = false;
if SpecAnalysis
    EpochLen = 1000;
    Overlap = 1000;
    ARC.Spectrumanalysis(ProjectPath ,'FreqBand', [1 40],'Conditions',{'REC1','REO'},'RedoifExist',true,'EpochLen' ,EpochLen ,'Overlap' ,Overlap,...
             'FigurePath','C:\Users\Elhamkhanom\Documents\My works\LongtermProject\ARCProject\Results\FFTData_NNRfitting',...
             'SavePath',SavePath);
end
%%

for S = 1:numel(SubIndex)
    % Read the subject's amplitude spectrum density, if already computed,
    % just loads them
    sub = SubIndex(S);
    display([SubjectData(sub).SubID]);
    temp = ARC.RES();
    if isempty(temp.loadRES(SavePath,SubjectData(sub)))
        error('RES class for subject not found, please run ARC.Spectrumanalysis first');
    else
        RESdata = temp.loadRES(SavePath,SubjectData(sub));
    end
    % extract REC and REO conditions
    Conds = RESdata.FindConditions({'REC1','REO'});
    FFTData = RESdata.FFTData(cell2mat(Conds));
    for i= 1:numel(FFTData)
        AllSubFFTData{S,i} = (mean(FFTData(i).Data));
    end
    Freq{S} = RESdata.Freq;
end

save(fullfile(SavePath,'REC_REO_Averagre_Specs'),'AllSubFFTData','Freq','SubjectData');

