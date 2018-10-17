function  model = ResParafac(RESdata,varargin)


opt = ParseArgs(varargin,...
    'FreqBand'      ,[1 15],...
    'Conditions'    ,[],...
    'Electrodes'    ,[],...
    'Corcondia'     , false,...
    'VarianceMode'  , 'spatial'...
    );

%% Select the conditions
if isempty(opt.Conditions)
    opt.Conditions = 'REC';
end
CondInd = RESdata.FindConditions(opt.Conditions);
CondInd = unique(cat(1,CondInd{:}));

FFTData = RESdata.FFTData(CondInd);
%% indicate the number of components using CorCondia
% Split half should be added here
InpData = cat(3,FFTData(:).Data);
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
end