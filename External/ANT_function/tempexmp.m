% Create some data: 

X = 1:12; 

Y = rand(1,12); 

% Generate a plot 

bar(X,Y); 

% Set the tick locations and remove the labels 

set(gca,'XTick',1:12,'XTickLabel','') 

% Define the labels 

lab = [{'January'};{'February'};{'March'};{'April'};{'May'};{'June'};...

           {'July'};{'August'};{'September'};{'October'};...

           {'November'};{'December'}];

% Estimate the location of the labels based on the position 

% of the xlabel 

hx = get(gca,'XLabel');  % Handle to xlabel 

set(hx,'Units','data'); 

pos = get(hx,'Position'); 

y = pos(2); 

% Place the new labels 

for i = 1:size(lab,1) 

    t(i) = text(X(i),y,lab(i,:)); 

end 

set(t,'Rotation',90,'HorizontalAlignment','right')  
