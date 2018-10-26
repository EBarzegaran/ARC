function Clusters = ClusterExtract(P,SupraTh,Stat,df,A)
% Extract clusters for cluster-based permutation test
% Inputs: 
    % P: is a Nx1 matrix consist of P values of statistical test
    % SupraTh: Threshold for extracting clusters
    % Stat: Nx1 test statistics from the test can be F or T value or any other
            % test stats
    % A: is a N x N matrix indicating the connectivity of the points (check ANOVAPermute for further help)
 
% OUTput:
    % Clusters: a cell array of cluster structures, they are sorted
            % descending according to their SStats
 
% Author: Elham Barzegaran
%% extract connected components
SigNod = find(P<SupraTh);
A2 = A(SigNod,SigNod);% subgraph consist of significant nodes
A2(1:length(A2)+1:end)=0;
G = graph(A2~=0);
bins = conncomp(G);

%% Prepare output
if ~isempty(bins)
    for b = 1:max(bins)
        Nodes{b} = SigNod(bins==b);
        SStat{b} = sum(Stat(SigNod(bins==b))); % different summary stats can be implemented here
    end
    Clusters = struct('Nodes',Nodes,'SStat',SStat);

    % sort clusters according to summaty stat
    [~,ord] = sort([Clusters.SStat],'descend');
    Clusters = Clusters(ord);
else
    Clusters = struct('Nodes',{[]},'SStat',{0});
end
end