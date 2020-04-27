function A = SourceNeighbors()
% Returns adjacency matrix which is Source x Source
% Author: Elham Barzegaran
%%
load('Source_CO_3005_nocel');
load('ds3000to400');

CO387 = CO(Sources,:);
Dist = squareform(pdist(CO387,'euclidean'));
Dist(1:length(Dist)+1:end) = 100;
A = double(Dist<7);% keep first order neighbors

end