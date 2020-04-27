function Source_visualization(data,direction,color,normal,hemis,Nodes)
% This function plots the sources on the surface
% INPUTS:
    % data: source data is a 5124*1, where 5124 is number of sources.
    % direction: the direction of headplot: it can has the values: 'up','down','left', 'right','all'
    % color: color of source: 'r', 'b', 'g', 'v', 'k' or a color map name like: 'jet'
    % normal: if the values should be normalized 1 or not 0;
    % hemis: it is either 'left' or 'righ', if no value is assigned , it will be considered as both
    
%% set default values
addpath('funcs');
if ~exist('direction','var'), direction ='up'; end
if ~exist('color','var'), color='r'; end
if ~exist('normal','var'), normal=1; end
if ~exist('hemis','var'), hemis  = [];end
if ~exist('Nodes','var'), Nodes  = [];end
ref_name = 'source_info';
%% plot units
%figure%
%('units','normalized','outerposition',[0 0 1 0.8])
if strcmp(direction,'all')
    subplot(2,3,3), sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); axis off;view(90,0);title('right') ;set(gca,'ylim',[-110 80]);set(gca,'Color',[0.8 0.8 0.8]);axis equal
    subplot(2,3,4), sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); set(gca,'Color',[0.8 0.8 0.8]);axis off;view(0,0);title('back');axis equal
    subplot(2,3,2), sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); set(gca,'Color',[0.8 0.8 0.8]);axis off;title('up');axis equal
    subplot(2,3,5), sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); set(gca,'Color',[0.8 0.8 0.8]);axis off;view(0,-90);title('down');axis equal
    subplot(2,3,1), sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); set(gca,'ylim',[-110 80]);set(gca,'Color',[0.8 0.8 0.8]);axis off;view(-90,0);title('left');axis equal
    subplot(2,3,6), sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); set(gca,'Color',[0.8 0.8 0.8]);axis off;view(180,0);title('front');axis equal
else
    sub_source_colored2(ref_name,data,color,normal,hemis,Nodes); axis off;

    switch direction
        case{'up'}

        case{'down'}
            view(90,-89)
            %lightangle(-180,180)
            
        case{'left'}
            view(-90,0)
        case{'right'}
            view(90,0)
        case{'front'}
            view(180,0)  
        case{'back'}
            view(0,0)
    end 
    set(gca,'Color',[0.8 0.8 0.8]);
end
axis equal
end