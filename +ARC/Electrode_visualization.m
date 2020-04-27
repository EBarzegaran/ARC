function Electrode_visulaization(Data,normali,color,eleclist)
% Data is a 64x1 vector including positive values(for now only positive)of power for each electrode
%normalize Data
% color: the desired color of map
if normali==1
    Data = (Data-min(Data(:)))/(max(Data(:))-min(Data(:)));
end
if ~exist('eleclist','var')
    eleclist = [];
end
Data = 1-Data;
%%
BadEl = [13 19];
V = importdata('elec_2D.txt');
V.data(BadEl,:)=[];V.textdata(BadEl)=[];
Data(BadEl)=[];
numEL=62;
%Dis1 = sqrt(sum((permute(repmat(V.data,[1 1 numEL]),[1 3 2])-permute(repmat(V.data,[1 1 numEL]),[3 1 2])).^2,3));
x = V.data(:,1);y = -V.data(:,2);%z=V.data(:,3);
tri=delaunay(x,y);
XX = ones(numel(x),3)*.9; 
switch color
    case 'r'
       XX(:,2:3)=repmat(Data,[1 2]);
    case 'g'
       XX(:,[1 3])=repmat(Data,[1 2]); 
    case 'b'
       XX(:,[1 2])=repmat(Data,[1 2]); 
    case 'm'
       XX(:,2)=Data; 
    case 'y'
       XX(:,3)=Data; 
    case 'c'
       XX(:,1)=Data; 
    case 'k'    
       XX(:,1:3)=repmat(Data,[1 3]);
    case 'hot'
        XX = 1-Data;
        colormap('hot');
    case {'hotcortex','coolhotcortex'}
       XX = 1-Data;
       colormap(jmaColors(color)); 
end
%patch('faces',tri,'vertices',[y x 0*x],'facevertexcdata',reshape(XX,numel(XX),1),'facecolor','interp','edgecolor','none');
patch('faces',tri,'vertices',[y x 0*x],'facevertexcdata',XX,'facecolor','interp','edgecolor','none');

hold on; 

if ~isempty(eleclist)
    Elecs = zeros(1,64);Elecs(eleclist)=1;Elecs = Elecs==1;
    Elecs(BadEl)=[];
    scatter(y(Elecs),x(Elecs),20,'b','linewidth',2,'Marker','o'); 
end
% load Cental_electrode.mat;EL = 1:64;EL(BadEl)=[];[E,I] = intersect (EL,Central_elec2);
% scatter(y(I),x(I),80,'r','filled'); 
%%
ph=linspace(0,2*pi,100)+pi/2; Coefx = 112; Coefy = 110;Coef = Coefx;
line(Coefx*cos(ph),Coefy*sin(ph)-10,'linewidth',1.5,'color','k','clipping','off');
n=3; line([cos(ph(n)) 1.15*cos(ph(1)) cos(ph(end-n+1))]*Coef,[sin(ph(n)) 1.15*sin(ph(1)) sin(ph(end-n+1))]*Coef-10,'linewidth',1.5,'color','k','clipping','off');
ph=linspace(pi-0.27*pi,pi+0.27*pi,20);
er=0.3;
ed=0.225;
line((-(1-ed)+er*cos(ph))*Coef,er*sin(ph)*Coef-10,'linewidth',1.5,'color','k','clipping','off');
line((+(1-ed)+er*cos(ph+pi))*Coef,er*sin(ph+pi)*Coef-10,'linewidth',1.5,'color','k','clipping','off');

xlim([-150 150]); ylim([-150 150]);
axis off;
axis equal;
end