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

StatResults = RmAnovaPermute(permute(ASDmat,[4 1 2 3]),A,100,.01,'mass');


%% Conduct post-hoc analysis
PHStatResults  = RmAnovaPostHoc(permute(ASDmat,[4 1 2 3]),StatResults);
Results.StatResults  = PHStatResults;

%% save the result
Results = []; 
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
    mt = 0;%floor(min(StatResults{i}.Uncorrected.F));
    Mt = ceil(max(StatResults{i}.Uncorrected.F));
    caxis([mt Mt]);
    
    CB = colorbar;
    set(get(CB,'title'),'string','F-Stats')
    CB.Ticks = mt:round((Mt-mt)/4):Mt;
    
    title(FacNames{i});
    set(S,'position',get(S,'position')+[-.05 -.05 .1 .1]);
    
end

set(Fhandler,'PaperPositionMode','manual');
set(Fhandler,'PaperPosition',[.25 .25 11 4]);
print(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.tif']),'-dtiff','-r300');
close;
%% plot the post-hoc results
levelNames{2} = opt.ModelNames;
levelNames{1} = {'ARC1','ARC2'};
for i = 1:numel(PHStatResults)
    FIG = figure;
    switch i
        case {1,2}
            Elecs = PHStatResults{i}.Clusters.Nodes;
            PHresults = PHStatResults{i}.Clusters.PostHoc;
            for ph = 1:numel(PHresults)
                Elec_F = zeros(64,1);
                Elec_F(Elecs) = PHresults(ph).F;
                SN = Elecs(PHresults(ph).P<(.01));
                SP = subplot(1,numel(PHresults),ph);
                [~,M1] = min(PHresults(ph).mean);%smaller level
                [~,M2] = max(PHresults(ph).mean);% larger level
                ARC.Electrode_visulaization(Elec_F,0,'hotcortex',SN);axis tight
                title([levelNames{PHresults(ph).factor}{PHresults(ph).levels(M2)} '>' levelNames{PHresults(ph).factor}{PHresults(ph).levels(M1)}]);
                caxis([0 max([PHresults.F])])
                %colorbar;
                set(SP,'position',get(SP,'position')+[-.005 -.005 .01 .01]);
            end
            set(FIG,'PaperPositionMode','manual');
            set(FIG,'PaperPosition',[1 1 4*numel(PHresults) 4]);
            
        case 3
            Elecs = PHStatResults{i}.Clusters.Nodes;
            PHresults = PHStatResults{i}.Clusters.PostHoc;
            % xaxis -> factor2
            % yaxis -> mean+SEM
            % two lines for factor 1
            ARC1M = [PHresults(1).mean PHresults(2).mean(2)];
            ARC1SEM = [PHresults(1).SEM PHresults(2).SEM(2)];

            ARC2M = [PHresults(4).mean PHresults(5).mean(2)];
            ARC2SEM = [PHresults(4).SEM PHresults(5).SEM(2)];
            
            FOIG = figure;
            hold on; p(1) = plot(1:3,ARC1M,'color',[0.1 .7 .1],'linewidth',2);
            errorbar(1:3,ARC1M,ARC1SEM,'.','color',[0.1 .7 .1],'linewidth',2);
            
            hold on; p(2) = plot([1:3]+.02,ARC2M,'color',[0.7 .1 .1],'linewidth',2);
            errorbar([1:3]+.02,ARC2M,ARC2SEM,'.','color',[0.7 .1 .1],'linewidth',2);
            % plot significant lines 
            
            set(gca,'xtick',1:3,'xticklabels',levelNames{2})
            xlim([.5 3.5])
            ylim([.5 4.5])
            ylabel('ASD')
            legend(p,levelNames{1})
            set(FIG,'PaperPosition',[1 1 7 3.5]);         
    end
    print(fullfile(opt.ResultsPath,'GroupLevel',['PostHoc_' FacNames{i} FileName '.tif']),'-dtiff','-r300');
    close all;
end
%% combine the figures

end



