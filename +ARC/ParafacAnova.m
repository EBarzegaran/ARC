function ParafacAnova(ModelPath,SubInfo,varargin)

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

StatResults = RmAnovaPermute(permute(ASDmat,[4 1 2 3]),A,100,.01,'mass');
%% plot the ANOVA results
FacNames = {'ARC','Condition','ARC X Condition'};
for i = 1:numel(StatResults)
    % Mark the significant clusters
    SC = [StatResults{i}.Clusters.Pvalue]<0.05;
    SN = [StatResults{i}.Clusters(SC).Nodes];
    %
    subplot(1,3,i),ARC.Electrode_visulaization(StatResults{i}.Uncorrected.F',0,'hotcortex',SN); axis tight
    colorbar;
    title(FacNames{i});
end

end



