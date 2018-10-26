function  PFresult = ResParafac(RESdata,varargin)

%Applied Parafac on RES data
%%
opt = ParseArgs(varargin,...
    'FreqBand'      ,[1 15],...
    'Conditions'    ,[],...
    'Electrodes'    ,[],...
    'Corcondia'     , false,...
    'VarianceMode'  , 'spatial',...
    'FixedFreqLoading', false,...
    'FixedModel'      ,[]...
    );
%% default values
if isempty(opt.Conditions)
    opt.Conditions = 'REC';
end
if isempty(opt.Electrodes)
    opt.Electrodes = 1:min(numel(RESdata.Elabels),64);
end

%% Select the conditions
CondInd = RESdata.FindConditions(opt.Conditions);
CondInd = unique(cat(2,CondInd{:}));

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
if opt.FixedFreqLoading
    MM = opt.FixedModel.Model;
    MM{1} = abs(rand(size(InpData,1),NumCom));
    MM{3} = abs(rand(size(InpData,3),NumCom));
end
if strcmpi(opt.VarianceMode,'temporal')
    if opt.FixedFreqLoading
        model = parafac(permute(InpData,[3 2 1]),NumCom,[0 0 0],[2 2 2],MM([3 2 1]),[0 1 0]);
    else
        model = parafac(permute(InpData,[3 2 1]),NumCom,[0 0 0],[2 2 2]);
    end
    model = model([3 2 1]);
elseif strcmpi(opt.VarianceMode,'spatial')
    if opt.FixedFreqLoading
        model = parafac(InpData,NumCom,[0 0 0],[2 2 2],MM,[0 1 0]);
    else
        model =  parafac(InpData,NumCom,[0 0 0],[2 2 2]);
    end
end

PFresult = ARC.PFModel(RESdata.SubjectInfo, model,{'Spatial','Frequency','Temporal'},opt.VarianceMode,...
    RESdata.Freq(FreqInd),RESdata.Elabels(opt.Electrodes),CondNames, CondLength);

PFresult =PFresult.OrganizeARCs(); % check the order of ARC1 and ARC2

end