function Fhandler = plotPFmodel(obj,isSave,Path)
% plot the loadings of a PF model
Fhandler=figure;
%% plot subject summary
h2=subplot(3,2,1);
TEXT = [obj.SubjectInfo.SubID ',' obj.SubjectInfo.Summaryinfo];
text(.25,.95,TEXT,'fontsize',12,'fontweight','bold'); axis off;
set(h2,'position',[0.33 0.85 0.1 0.1]);
colors = [0 1 0;1 0 0; 0 0 1; 0 1 1 ; 0 0 0; 1 1 0; 1 0 1] ;
Cleg = {'ARC1','ARC2','ARC3'};
FS = 10;
%% plot frequency loading
h2=subplot(3,3,1+3); hold on;

Comp =  obj.SubjectInfo.Comp;Comp = Comp(Comp>0);
NumCom = numel(Comp);
for C = 1:NumCom
    plot(obj.Freq,obj.Model{strcmpi(obj.Modes,'frequency')}(:,Comp(C)),'Color',colors(C,:) ,'LineWidth',1.5);
end 
title('Frequency Loadings','FontSize',FS,'Fontweight','bold');
xlabh1=xlabel('Frequency (Hz)','FontSize',FS,'Fontweight','bold');
xlim([6 14]);
set(gca,'FontSize',FS-2); 
L=legend(Cleg{1:NumCom});
set(L,'position',[.1 .75 .1 .1])
set(h2,'position',[0.1 0.50 0.30 0.35]);
%set(xlabh,'Position',get(xlabh,'Position') - [0 .005 0]);

%% Plot mean temporal loadings

h2= subplot(3,3,2+3); 
[TMeans, TSEM] = obj.PFTempMeans;
bar_h = bar(TMeans(:,Comp)); hold on;
if numel(bar_h)==1
   errorbar(bar_h.XData+bar_h.XOffset,TMeans(Comp),TSEM(Comp),'.','color','k');
   bar_h.FaceColor = 'flat';
   for c = 1:NumCom
       bar_h.CData(c,:) = colors(c,:);
   end
else
    for c = 1:NumCom
        errorbar(bar_h(c).XData+bar_h(c).XOffset,TMeans(:,Comp(c)),TSEM(:,Comp(c)),'.','color','k');
        bar_h(c).FaceColor = colors(c,:);
    end
end
title('Temporal Loadings','FontSize',FS,'Fontweight','bold');
set(gca,'FontSize',FS-2);
xlabh2= xlabel('Conditions','FontSize',FS,'Fontweight','bold');

set(gca,'XTick',1:numel(obj.Conditions),'XTickLabel',obj.Conditions,'FontSize',FS-2);
rotateXLabels( gca(), 45 );box off
set(h2,'position',[0.55 0.50 0.33 0.35]);
if NumCom==1
    set(h2,'position',[0.55 0.50 0.23 0.35]);
end
% P2 = get(xlabh2,'Position');
% P1= get(xlabh1,'Position');
% set(xlabh1,'Position',[P1(1) P2(2) P1(3)]);
%% Plot spatial loadings
colors = {'g','r','b','c'};
load('ElectrodeLabels.mat');
for C = 1:NumCom
    modelallelec = zeros(64,1);
    % find the electrodes
    Elecs = cellfun(@(x) find(strcmpi(Elabel,x)),obj.Elabel,'uni',false);
    Elecs = cat(1,Elecs{:});
    modelallelec(Elecs) = obj.Model{strcmpi(obj.Modes,'spatial')}(:,Comp(C));
    h2=subplot(3,4,C+8); ARC.Electrode_visulaization(modelallelec,1,colors{C});title(Cleg{C},'FontSize',FS,'Fontweight','bold');
    set(h2,'position',[0.05+(C-1)*.23 0.02 0.20 0.30]);
end

%% Save the figure
set(Fhandler,'PaperPositionMode','manual');
set(Fhandler,'PaperPosition',[.25 .25 9 7 ]);
if isSave
    print(fullfile(Path,[obj.SaveName '.tif']),'-dtiff','-r300');
end
end