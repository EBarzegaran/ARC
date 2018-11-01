function Clusters = ClusterExtract(P,SupraTh,Stat,BaseF,A,StatType)
% Extract clusters for cluster-based permutation test
% Inputs: 
    % P: is a Nx1 matrix consist of P values of statistical test
    % SupraTh: P-value Threshold for extracting clusters
    % Stat: Nx1 test statistics from the test can be F or T value or any other
            % test stats
    %BaseF: is the Fstat calculated as the f-inverse of the SupraTh
    % A: is a N x N matrix indicating the connectivity of the points (check ANOVAPermute for further help)
    %StatType: ['mass']/'size'/'height'/'TFCE': check: Pernet, C. R., et al. "Cluster-based computational methods for mass univariate analyses of event-related brain potentials/fields: A simulation study." Journal of Neuroscience Methods 250 (2015): 85-93.
 
% OUTput:
    % Clusters: a cell array of cluster structures, they are sorted
            % descending according to their SStats
 
% Author: Elham Barzegaran

%%
if ~exist('StatType','var')
    StatType = 'mass';
end
%%
if strcmp(StatType,'TFCE')
    H = 1;E = .5;
    Mf = max(Stat);mf = min(Stat);
    steps = mf:(Mf-mf)/1000:Mf;
    TFCE = zeros(size(Stat));
    for s = 1:numel(steps)
        h = steps(s);% height
        % calculating cluster sizes
        SigNod = find(Stat>=h);
        A2 = A(SigNod,SigNod);% subgraph consist of significant nodes
        A2(1:length(A2)+1:end)=0;% diagonal to zero
        G = graph(A2~=0);
        bins = conncomp(G);% connected components
        CS = arrayfun(@(x) sum(bins==x),1:max(bins));% cluster sizes
        CSN = zeros(size(bins));
        for b = 1: max(bins), CSN(bins==b) = CS(b);end
        e = zeros(size(Stat));
        e(SigNod) = CSN;
        TFCE = TFCE + (e.^E)*(h^H);
    end
    TFCE = TFCE/1000;
    error('TFCE is under development. Please try another cluster statistics')
else
    % extract connected components
    SigNod = find(P<SupraTh);
    A2 = A(SigNod,SigNod);% subgraph consist of significant nodes
    A2(1:length(A2)+1:end)=0;
    G = graph(A2~=0);
    bins = conncomp(G);

    % Prepare output
    if ~isempty(bins)
        for b = 1:max(bins)
            Nodes{b} = SigNod(bins==b);
            switch StatType
                case 'mass'
                    SStat{b} = sum(Stat(SigNod(bins==b))-BaseF); % different summary stats can be implemented here
                case 'size'
                    SStat{b} = numel(Nodes{b});
                case 'height'
                    SStat{b} = max(Stat(SigNod(bins==b)));
                otherwise
                    error('The stat type is not found. Please check ClusterExtract help.');
            end
        end
        Clusters = struct('Nodes',Nodes,'SStat',SStat);

        % sort clusters according to summaty stat
        [~,ord] = sort([Clusters.SStat],'descend');
        Clusters = Clusters(ord);
    else
        Clusters = struct('Nodes',{[]},'SStat',{0});
    end
end
end