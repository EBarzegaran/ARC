function Results = ParafacAnova(ModelPath,SubInfo,varargin)

%% Parse input and assign default values
opt = ParseArgs(varargin,...
    'FileNames'     ,[],...
    'ModelNames'    ,[],...
    'ResultsPath'   ,[],...
    'VarianceMode'  ,'temporal' ...
    );
%% Load PARAFAC models
for sub = 1:numel(SubInfo)
    for ml = 1:numel(opt.FileNames)
        ModelTemp = ARC.PFModel();
        ModelTemp = ModelTemp.loadPFModel(ModelPath,...
            SubInfo(sub),opt.FileNames{ml},opt.VarianceMode);
        ModelAll{sub,ml} = ModelTemp;
    end
end
% get loadings 
FreqLoad = cellfun(@(x) x.GetLoading([1 2],'frequency'),ModelAll,'uni',false);
SpatLoad = cellfun(@(x) x.GetLoading([1 2],'spatial'),ModelAll,'uni',false);
TempLoad = cellfun(@(x) x.GetLoading([1 2],'temporal'),ModelAll,'uni',false);
% calculate ASDs 
ASD = cellfun(@(x,y,z) repmat(max(x).*mean(z),[size(y,1) 1]).*y, FreqLoad,SpatLoad,TempLoad,'uni',false);

%% ANOVA, main effect and interactions
ASDmat = arrayfun(@(x) cat(3,ASD(x,:)), 1:size(ASD,1),'uni',false);
ASDmat = cellfun(@(x) cat(3,x{:}),ASDmat,'uni',false );ASDmat = cat(4,ASDmat{:});
A = ElectrodeNeighbors();

StatResults = RmAnovaPermute(permute(ASDmat,[4 1 2 3]),A,1000,.01,'mass');

%% save the result
Results = []; 
Results.StatResults = StatResults;
Results.SubjectInfo = SubInfo;
Results.Conditions = opt.ModelNames;
Results.VarianceMode = opt.VarianceMode;
if ~exist(fullfile(opt.ResultsPath,'GroupLevel'),'dir')
    mkdir(fullfile(opt.ResultsPath,'GroupLevel'));
end
FileName = ['GroupLevelANOVA_' [Results.Conditions{:}] '_Var' Results.VarianceMode '_' num2str(numel(Results.SubjectInfo)) 'Subs'];
save(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.mat']),'Results');

%% plot the ANOVA results
FacNames = {'ARC','Condition','ARC X Condition'};
Fhandler = figure;
for i = 1:numel(StatResults)
    % Mark the significant clusters
    SC = [StatResults{i}.Clusters.Pvalue]<0.05;
    SN = [StatResults{i}.Clusters(SC).Nodes];
    %
    S = subplot(1,3,i);
    ARC.Electrode_visulaization(StatResults{i}.Uncorrected.F',0,'hotcortex',SN); axis tight
    CB = colorbar;
    CBP = CB.Position; CBP(3) = .03; CB.Position = CBP;
    set(get(CB,'title'),'string','F-Stats')
    title(FacNames{i});
end

set(Fhandler,'PaperPositionMode','manual');
set(Fhandler,'PaperPosition',[.25 .25 12 3.5]);
print(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.tif']),'-dtiff','-r300');
end



