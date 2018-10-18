function  PFresult = ResParafac(RESdata,varargin)


opt = ParseArgs(varargin,...
    'FreqBand'      ,[1 15],...
    'Conditions'    ,[],...
    'Electrodes'    ,[],...
    'Corcondia'     , false,...
    'VarianceMode'  , 'spatial'...
    );
%% default values
if isempty(opt.Conditions)
    opt.Conditions = 'REC';
end
if isempty(opt.Electrodes)
    opt.Electrodes = 1:numel(RESdata.Elabels);
end

%% Select the conditions
CondInd = RESdata.FindConditions(opt.Conditions);
CondInd = unique(cat(1,CondInd{:}));

FFTData = RESdata.FFTData(CondInd);
InpData = cat(3,FFTData(:).Data);

FreqInd = (RESdata.Freq>=opt.FreqBand(1) & RESdata.Freq<=opt.FreqBand(2));
InpData = InpData(opt.Electrodes,FreqInd,:);
%% Select frequencyband and Electrodes
CondNames = RESdata.GetCondNames(CondInd);
CondLength = [RESdata.FFTData(CondInd).Epochs];

%% indicate the number of components using CorCondia
% Split half should be added here
if opt.Corcondia
    [~,Corco,~] = pftest(3,InpData,4,[0 0 0],[2 2 2]);%
    NumCom = find(mean(Corco,2)>80,1,'last');
else
    NumCom = RESdata.SubjectInfo.Compnum;
end

%% Apply Parafac

if strcmpi(opt.VarianceMode,'temporal')
    model = parafac(permute(InpData,[3 2 1]),NumCom,[0 0 0],[2 2 2]);
    model = model([3 2 1]);
elseif strcmpi(opt.VarianceMode,'spatial')
    model =  parafac(InpData,NumCom,[0 0 0],[2 2 2]);
end

PFresult = ARC.PFModel(RESdata.SubjectInfo, model,{'Spatial','Frequency','Temporal'},...
    RESdata.Freq(FreqInd),RESdata.Elabels(opt.Electrodes),CondNames, CondLength);
end