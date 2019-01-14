function Spectrumanalysis(ProjectPath,varargin)
% This function reads .cnet files, calculate spectrums and save them as RES class
% varargin (optional):
    % Conditions: Note that this is only for visualization (Spectrum plots), RES files contain
                  % all the conditions

opt = ParseArgs(varargin,...
    'EpochLen'   ,2500,...
    'Overlap'    ,1250,...
    'FreqBand'  ,[1 15],...
    'Space'      ,'Electrode',...
    'Conditions' ,[],...
    'Subjectinfo'   ,[],...
    'SubjectSelect' ,[],...
    'FigurePath'   ,[],...
    'RedoifExist'   ,false,...
    'PlotResults'   ,true,...
    'SavePath'      ,fullfile(ProjectPath ,'FFTData')... 
        );
    
EpLen = opt.EpochLen; % Epoch length
MovWin = opt.Overlap;%Ovelpa of the moving window

if isempty(opt.Subjectinfo)
    load('+ARC\Private\Subjectinfo.mat');
    opt.Subjectinfo = SubjectData;
end

if isempty(opt.SubjectSelect)% select the first recordings of all subjects
    opt.SubjectSelect = find([SubjectData(:).Longitude]==0);
end

%% Read files, calculate spectrums and save them as RES class
if ~exist(fullfile(ProjectPath ,'FFTData'),'dir')
    mkdir(fullfile(ProjectPath ,'FFTData'));
end
for S = 1:numel(opt.SubjectSelect)
    sub = opt.SubjectSelect(S);
    display([SubjectData(sub).SubID]);
    temp = ARC.RES();
    if isempty(temp.loadRES(opt.SavePath,SubjectData(sub))) || opt.RedoifExist
        [RESdata] = RawEEGtoRES(ProjectPath,SubjectData(sub),opt.EpochLen ,opt.Overlap ,opt.FreqBand(2),opt.FreqBand(1));
        RESdata.saveRES(opt.SavePath);
    else
        RESdata = temp.loadRES(opt.SavePath,SubjectData(sub));
    end
    if opt.PlotResults
        if ~exist(fullfile(opt.FigurePath,'Specrum'),'dir')
            mkdir(fullfile(opt.FigurePath,'Specrum'));
        end
        RESdata.PlotSpectrum('SavePath',fullfile(opt.FigurePath,'Specrum'),'Conditions',opt.Conditions,'FreqBand',opt.FreqBand);
        close;
    end
end

end