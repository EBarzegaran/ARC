function A = ElectrodeNeighbors()
% Returns adjacency matrix which is Electrode x Electrodes
% Author: Elham Barzegaran
%%
V = importdata('elec_2D.txt');
x = V.data(:,1);y = -V.data(:,2);%z=V.data(:,3);
tri=delaunay(x,y);

%% 
A = zeros(numel(x));
for x = 1:size(tri,1)
    A(tri(x,:),tri(x,:)) = 1;
end
end