function StatResults = ParafacAnova(ModelPath,SubInfo,varargin)

%% Parse input and assign default values
opt = ParseArgs(varargin,...
    'FileNames'     ,[],...
    'ModelNames'    ,[],...
    'ResultsPath'   ,[],...
    'VarianceMode'  ,'temporal', ...
    'Space'         ,'Electrode',...
    'redoAnalysis'  , true,...
    'redoANOAVfigs' ,true,...
    'PermNum'       ,1000 ...
    );
%% Load PARAFAC models
for sub = 1:numel(SubInfo)
    for ml = 1:numel(opt.FileNames)
        ModelTemp = ARC.PFModel();
        ModelTemp = ModelTemp.loadPFModel(ModelPath,...
            SubInfo(sub),opt.FileNames{ml},opt.VarianceMode,opt.Space);
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
if strcmp(opt.Space,'Electrode')
    A = ElectrodeNeighbors();
elseif strcmp(opt.Space,'Source')
    A = SourceNeighbors();
else
    error('The space can be either electrode or source');
end

FileName = ['GroupLevelANOVA_' [opt.ModelNames{:}] '_Var' opt.VarianceMode '_' num2str(numel(SubInfo)) 'Subs_' opt.Space];
if opt.redoAnalysis || ~exist(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.mat']),'file')
    StatResults = RmAnovaPermute(permute(ASDmat,[4 1 2 3]),A,opt.PermNum,.01,'mass');

    if strcmp(opt.Space,'Source')
        for i = 1:numel(StatResults)
            for j = 1:numel(StatResults{i}.Clusters)
                [brod_label, lobe_label, hem_label, aal_label]= ParcellateSource(StatResults{i}.Clusters(j).Nodes);
                StatResults{i}.Clusters(j).Labels = [num2cell(brod_label) lobe_label  hem_label aal_label];
            end
        end
    end

    %% Conduct post-hoc analysis
    PHStatResults  = RmAnovaPostHoc(permute(ASDmat,[4 1 2 3]),StatResults);
    Results.StatResults  = PHStatResults;

    %% save the result
    %Results = []; 
    Results.SubjectInfo = SubInfo;
    Results.Conditions = opt.ModelNames;
    Results.VarianceMode = opt.VarianceMode;
    Results.Space = opt.Space;
    if ~exist(fullfile(opt.ResultsPath,'GroupLevel'),'dir')
        mkdir(fullfile(opt.ResultsPath,'GroupLevel'));
    end
    %FileName = ['GroupLevelANOVA_' [Results.Conditions{:}] '_Var' Results.VarianceMode '_' num2str(numel(Results.SubjectInfo)) 'Subs_' Results.Space];
    save(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.mat']),'Results');
else
    load(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.mat']),'Results');
    PHStatResults = Results.StatResults;
end

%% plot the ANOVA results
FacNames = {'ARC','Condition','ARC X Condition'};
if opt.redoANOAVfigs
    if strcmp(opt.Space,'Source')
        Orient = {'left','back','right'};
    else
        Orient = 'topo';
    end
    noClust = false;
    Fhandler = figure;
    
    set(Fhandler,'PaperPositionMode','manual');
    set(Fhandler,'PaperPosition',[.25 .25 7 2*numel(Orient)]);
    set(Fhandler,'Unit','inch','Position',[0 0 7 2*numel(Orient)]);
    for Or = 1:numel(Orient)
        for i = 1:numel(PHStatResults)
            % Mark the significant clusters
            SC = [PHStatResults{i}.Clusters.Pvalue]<0.05;
            SN = [PHStatResults{i}.Clusters(SC).Nodes];
            if strcmp(opt.Space,'Electrode')
                S = subplot(1,3,i);
                if noClust
                    ARC.Electrode_visulaization(PHStatResults{i}.Uncorrected.F',0,'hotcortex',SN); axis tight
                else
                    ARC.Electrode_visulaization(PHStatResults{i}.Uncorrected.F',0,'hotcortex',[]); axis tight
                end
                mt = 0;%floor(min(PHStatResults{i}.Uncorrected.F));
                 Mt = ceil(max(PHStatResults{i}.Uncorrected.F));
                 axis tight; axis vis3d
                 SP = get(S,'position');
                title(FacNames{i});
                CB = colorbar;
                set(get(CB,'title'),'string','F-Stats')
                set(S,'position',SP+[-.05 -.05 -.05 -.05]);
                
            elseif strcmp(opt.Space,'Source')
                S = subplot(3,3,(i-1)*3+Or);
                load('MaptoSurface_400to5000.mat');
                Th =0.3;
                m1=PHStatResults{i}.Uncorrected.F*D;
                %m1=(m1-min(m1(:)))/(max(m1(:))-min(m1(:)));
                [~,MapInd] = (max(D'));
                if noClust
                    ARC.Source_visualization(m1,Orient{Or},'hotcortex',0,[],[]);%MapInd(SN));
                else 
                    ARC.Source_visualization(m1,Orient{Or},'hotcortex',0,[],MapInd(SN));
                end
            
             mt = 0;%floor(min(PHStatResults{i}.Uncorrected.F));
             Mt = ceil(max(PHStatResults{i}.Uncorrected.F));
             axis tight; axis vis3d
             SP = get(S,'position');
             if Or ==3
                CB = colorbar;
                set(get(CB,'title'),'string','F-Stats')
                %CB.YAxisLocation = 'left';
                set(CB,'position',get(CB,'position')+[.06 0.01 -.00 -.08]);
             end
            if Or==1 
                if i<3
                    axes('Position',[0.05 0.72-(.3*(i-1)) 0.05 0.15]); axis off
                else
                    axes('Position',[0.05 0.1 0.05 0.15]); axis off
                end
                text(0,0,FacNames{i},'fontsize',10,'rotation',90);
                %axis on;
            end
            if Or~=2
                set(S,'position',SP+[-.07 -.07 .05 .05]);
            else
                set(S,'position',SP+[-.04 .0 .0 .0]);
            end
            
            end
            
        end
    end

    %
    if noClust
        print(fullfile(opt.ResultsPath,'GroupLevel',[FileName 'NoCluster.tif']),'-dtiff','-r300');
    else
        print(fullfile(opt.ResultsPath,'GroupLevel',[FileName '.tif']),'-dtiff','-r300');
    end
    close;
end
%% plot the post-hoc results, we plot the post-hoc results separately for each cluster
levelNames{2} = opt.ModelNames;
levelNames{1} = {'ARC1','ARC2'};
Cols = [0.1 .7 .1; 0.7 .1 .1];    
for i = 1:numel(PHStatResults)

    switch i
        case {1,2} % main effects
            for cl = 1:numel(PHStatResults{i}.Clusters)
                if ~isempty(PHStatResults{i}.Clusters(cl).PostHoc)

                    FIG = figure;
                    Nodes = PHStatResults{i}.Clusters(cl).Nodes;
                    PHresults = PHStatResults{i}.Clusters(cl).PostHoc;
                    for ph = 1:numel(PHresults)
                        Nodes_F = zeros(size(ASD{1,1},1),1);
                        Nodes_F(Nodes) = PHresults(ph).F;
                        SN = Nodes(PHresults(ph).P<(.01));
                        SP = subplot(1,numel(PHresults),ph);
                        [~,M1] = min(PHresults(ph).mean);%smaller level
                        [~,M2] = max(PHresults(ph).mean);% larger level
                        if strcmp(opt.Space,'Electrode')
                            ARC.Electrode_visulaization(Nodes_F,0,'hotcortex',SN);axis tight
                        elseif strcmp(opt.Space,'Source')
                            load('MaptoSurface_400to5000.mat');
                            [~,MapInd] = (max(D'));
                            m1=Nodes_F'*D;
                            ARC.Source_visualization(m1,'back','hotcortex',0,[],MapInd(SN));
                        end
                        title([levelNames{PHresults(ph).factor}{PHresults(ph).levels(M2)} '>' levelNames{PHresults(ph).factor}{PHresults(ph).levels(M1)}]);
                        caxis([0 max([PHresults.F])])
                        colorbar;
                        if strcmp(opt.Space,'Electrode')
                            set(SP,'position',get(SP,'position')+[-.005 -.005 .01 .01]);
                        end
                    end
                    set(FIG,'PaperPositionMode','manual');
                    set(FIG,'PaperPosition',[1 1 4*numel(PHresults) 4]);
                    print(fullfile(opt.ResultsPath,'GroupLevel',['PostHoc_' FacNames{i} FileName '_Cluster' num2str(cl) '.tif']),'-dtiff','-r300');
                    close all;
                end
            end


        case 3 % plot interactions
            SigClust = find(cellfun(@(x) ~isempty(x),{PHStatResults{i}.Clusters.PostHoc}));
            Clusnum = numel(SigClust);
            FIG = figure;
            for c = 1:Clusnum%numel(PHStatResults{i}.Clusters)
                cl = SigClust(c);
                if ~isempty(PHStatResults{i}.Clusters(cl).PostHoc)
                    subplot(1,Clusnum,c)
                    % which hemisphere?
                    hemi = PHStatResults{i}.Clusters(cl).Labels{1,5};
                    if strcmp(hemi,'L')
                        title('Left cluster')
                    elseif strcmp(hemi,'R')
                        title('Right cluster')
                    end
                        
                    
                    Elecs = PHStatResults{i}.Clusters(cl).Nodes;
                    PHresults = PHStatResults{i}.Clusters(cl).PostHoc;

                    % two lines for factor 1
                    if numel(levelNames{2})>2
                        ARC1M = [PHresults(1).mean PHresults(2).mean(2)];
                        ARC1SEM = [PHresults(1).SEM PHresults(2).SEM(2)];
                        ARC2M = [PHresults(4).mean PHresults(5).mean(2)];
                        ARC2SEM = [PHresults(4).SEM PHresults(5).SEM(2)];
                    else
                        ARC1M = [PHresults(1).mean];
                        ARC1SEM = [PHresults(1).SEM];
                        ARC2M = [PHresults(2).mean];
                        ARC2SEM = [PHresults(2).SEM];
                    end
                    hold on; p(1) = plot(1:numel(levelNames{2}),ARC1M,'color',Cols(1,:),'linewidth',2);
                    errorbar(1:numel(levelNames{2}),ARC1M,ARC1SEM,'.','color',Cols(1,:),'linewidth',2);

                    hold on; p(2) = plot([1:numel(levelNames{2})]+.02,ARC2M,'color',Cols(2,:),'linewidth',2);
                    errorbar([1:numel(levelNames{2})]+.02,ARC2M,ARC2SEM,'.','color',Cols(2,:),'linewidth',2);

                    % plot significant lines 
                    M = max(max(ARC1M+ARC1SEM),max(ARC2M+ARC2SEM));
                    m = min(min(ARC1M-ARC1SEM),min(ARC2M-ARC2SEM));
                    R = M-m;
                    for comp = 1:2*nchoosek(numel(levelNames{2}),2)
                        if mean(PHStatResults{3}.Clusters(cl).PostHoc(comp).P)<.1
                            ARCnum = PHStatResults{3}.Clusters(cl).PostHoc(comp).levels1;
                            if ARCnum ==2
                                y = M + (R*.05)*(5+ ((comp-1-(numel(levelNames{2}))))*1.7);
                                y2 = y+ R*.035;
                                y3 = y- R*.03;
                            else
                                y = m - (R*.05)*(4 - ((comp-1))*1.7);
                                y2 = y - R*.03;
                                y3 = y + R*.03;
                            end

                            x = PHStatResults{3}.Clusters(cl).PostHoc(comp).levels2;
                            line(x,[y y],'LineStyle','-.','color',Cols(ARCnum,:),'linewidth',1.2);
                            %text(mean(x),y2,'*')
                            text(mean(x)-.2,y2,['P = ' num2str(round(mean(PHStatResults{3}.Clusters(cl).PostHoc(comp).P),3))],'fontsize',8)
                            line([x(1) x(1)],[y y3],'LineStyle','-.','color',Cols(ARCnum,:),'linewidth',1.2);
                            line([x(2) x(2)],[y y3],'LineStyle','-.','color',Cols(ARCnum,:),'linewidth',1.2);
                        end
                    end

                    for comp = 2*nchoosek(numel(levelNames{2}),2)+1:2*nchoosek(numel(levelNames{2}),2)+numel(levelNames{2})
                        if mean(PHStatResults{3}.Clusters(cl).PostHoc(comp).P)<.1
                            x = PHStatResults{3}.Clusters(cl).PostHoc(comp).levels1;
                            line([x x]-.1,[ARC1M(x) ARC2M(x)],'linewidth',1.2,'LineStyle','-.','color','k');
                            line([x-.1 x-.07],[ARC1M(x) ARC1M(x)],'linewidth',1.2,'LineStyle','-.','color','k');
                            line([x-.1 x-.07],[ARC2M(x) ARC2M(x)],'linewidth',1.2,'LineStyle','-.','color','k');
                            text(x-.15,(ARC1M(x)+ARC2M(x))/2.14,['P = ' num2str(round(mean(PHStatResults{3}.Clusters(cl).PostHoc(comp).P),2))],'fontsize',8,'rotation',90)
                        end
                    end

                    set(gca,'xtick',1:3,'xticklabels',levelNames{2})
                    xlim([.5 numel(levelNames{2})+.5])
                    ylim([m*.6 M*1.33])
                    ylabel('ASD')
                    L = legend(p,levelNames{1});
                    set(L,'position',get(L,'position')+[0 0.04 0 0]);
                    
                end
            end
            set(FIG,'unit','inch','PaperPosition',[1 1 5*Clusnum 5],'Position',[1 1 5*Clusnum 5]); 
            print(fullfile(opt.ResultsPath,'GroupLevel',['PostHoc_' FacNames{i} FileName '.tif']),'-dtiff','-r300');
            close all;
    end

end

end



