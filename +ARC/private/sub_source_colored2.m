function sub_source_colored2(ref_name,source_data,Color,normal,hemis,Nodes)

load(ref_name);
%g3.vertices(:,1) = -g3.vertices(:,1);

D2 = prepare_refine(source_data,g,ind1,val1,size(g.vertices,1),size(g2.vertices,1));
D3 = prepare_refine(D2,g2,ind2,val2,size(g2.vertices,1),size(g3.vertices,1));
data = D3;
if ~isempty(Nodes)
    scatter3(g3.vertices(Nodes,1),g3.vertices(Nodes,2),g3.vertices(Nodes,3),10,'g','filled');
    hold on;
end

%% prepare values
% Red color
PD = data(data>=0);PI = find(data>=0);
if Color == 'r'
    if normal ==1
        PD = ceil((PD-min(PD))/(max(PD)-min(PD))*31);
    else
        PD = ceil(PD*31);
    end
    if isnan(mean(PD)),PD(isnan(PD)) = .0;PD = ceil(PD*31);end

    PD = PD+1;
    cmap = zeros(32,3); % color map postive
    cmap(:,1) = 0.901:.1/32:1;
    cmap(:,2) = 0.891:-.9/32:.0;
    cmap(:,3) = 0.891:-.9/32:.0;
    XX(PI,:) = cmap(PD,:);
elseif Color == 'b'
        % Blue color
    if normal ==1
        PD = ceil((PD-min(PD))/(max(PD)-min(PD))*31);
    else
    PD = ceil(PD*31);
    end
    if isnan(mean(PD)),PD(isnan(PD)) = .0;PD = ceil(PD*31);end

    PD = PD+1;
    cman = zeros(32,3); % color map negative
    cman(:,1) = 0.891:-.9/32:.0;
    cman(:,2) = 0.891:-.9/32:.0;
    cman(:,3) = 0.901:.1/32:1;
    XX(PI,:) = cman(PD,:);
elseif Color == 'g'
% Green color
    if normal ==1
        PD = ceil((PD-min(PD))/(max(PD)-min(PD))*31);
    else
        PD = ceil(PD*31);
    end
    if isnan(mean(PD)),PD(isnan(PD)) = .0;PD = ceil(PD*31);end

    PD = PD+1;
    cman = zeros(32,3); % color map negative
    cman(:,1) = 0.891:-.9/32:.0;
    cman(:,2) = 0.802:.1/32:.901;
    cman(:,3) = 0.891:-.9/32:.0;
    XX(PI,:) = cman(PD,:);
elseif Color == 'v'
% Green color
    if normal ==1
        PD = ceil((PD-min(PD))/(max(PD)-min(PD))*31);
    else
        PD = ceil(PD*31);
    end
    if isnan(mean(PD)),PD(isnan(PD)) = .0;PD = ceil(PD*31);end

    PD = PD+1;
    cman = zeros(32,3); % color map negative
    cman(:,1) = 0.802:.1/32:.901;
    cman(:,2) = 0.891:-.9/32:.0;
    cman(:,3) = 0.802:.1/32:.901;
    XX(PI,:) = cman(PD,:);
elseif Color == 'k'
% Green color
    PD = ceil((PD-min(PD))/(max(PD)-min(PD))*31);
    if isnan(mean(PD)),PD(isnan(PD)) = .0;PD = ceil(PD*31);end

    PD = PD+1;
    cman = zeros(32,3); % color map negative
    cman(:,1) = 0.891:-.9/32:.0;
    cman(:,2) = 0.802:.1/32:.901;
    cman(:,3) = 0.802:.1/32:.901;
    XX(PI,:) = cman(PD,:);
else
    if normal ==1
        %PD = ceil((PD-min(PD))/(max(PD)-min(PD))*31);
        PD = ((PD-min(PD))/(max(PD)-min(PD)));
    else
        %PD = ceil(PD*31);
        PD = ((PD-min(PD))/(max(PD)-min(PD)))*max(source_data);
    end
    %if isnan(mean(PD)),PD(isnan(PD)) = .0;PD = ceil(PD*31);end
    %PD = PD+1;
    
    XX = PD';
end

if strcmp(hemis,'left')||strcmp(hemis,'right')
    [ind g3] = selec_hemisphere(g3,hemis);
    XX = XX(ind,:);
end
%%
XX(isnan(XX))=0;

patch('faces',g3.faces,'vertices',g3.vertices,'edgecolor','none','facecolor','interp','facevertexcdata',XX,'FaceAlpha',0.85,...
    'AmbientStrength',0.45,'SpecularStrength',0.1,'DiffuseStrength', 0.45,'SpecularExponent',10,'FaceLighting','gouraud');axis off;%set(gca,'Color',[.1 .1 1]);
if sum(strcmpi(Color,{'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'}))
    colormap(Color);
elseif sum(strcmpi(Color,{'arizona','Air force','USC','Cal','Nebraska','italy','hotcortex','coolhotcortex','pval','phasecolor','coolhot'}))
    colormap(jmaColors(Color));
end

shading interp
lightangle(90,10)
lightangle(-90,10)
%lightangle(90,-50)
end
