function Clusters = ANOVAPermute(Data,A,PermNum,SupraTh)
    % This function runs a cluster based permuation test for within-subject
    % repeated measures ANOVA or paired TTEST
    % INPUT:
        % Data: Subjects X Electrodes(Time or Source) X Factor1 X Factor2  is
                % the input data, if you want to do TTest then the size of Factor2=1;
        % A: is the adjacency matrix which should be Electrode X Electrode
                % or Source X Source or Time X Time matrix, indicting the
                % neighboring electrodes or sources. If the first dimension
                % of Data is time, then A = diag(ones(Time,1),-1)+ diag(ones(Time,1),1)
                
        % permnum: is the number of permutations, minimum 1000 permutations
                % is recommended
                
   % Author: Elham Barzegaran, 10/2018
   %% default values
   if ~exist('PermNum','var') || isempty(PermNum)
       PermNum = 1000;
   end
   
   if ~exist('SupraTh','var') || isempty(SupraTh)
       SupraTh = .05;
   end
   
   %% Prepare the data structure and teh grouping labels
    Y = squeeze(Data(:,1,:,:));
    SubID = repmat((1:size(Y,1))',[1 size(Y,2) size(Y,3)]);
    Fac1 = permute(repmat((1:size(Y,2))',[1 size(Y,1) size(Y,3)]),[2 1 3]);% ARC
    Fac2 = permute(repmat((1:size(Y,3))',[1 size(Y,2) size(Y,1)]),[3 2 1]);% cond

    %% cluster based-permuation test ANOVA 
    Factors = {'fac1','fac2','fac1 x fac2'};
    for perm = 1: PermNum+1
        if mod(perm,20)==0,disp(num2str(perm));end
        if perm==1 % unpermuted data
            % Calculate ANOVA
            for El = 1:size(Data,2)
                Y = squeeze(Data(:,El,:,:));
                stats{El} = rm_anova2(Y(:),SubID(:),Fac1(:),Fac2(:),Factors(1:2));
            end
            % extract cluster statistics
            for fac = 1:3
                P = cellfun(@(x) x(strcmpi(stats{1}(:,1),Factors{fac}),strcmpi(stats{1}(1,:),'P')),stats);P = [P{:}];
                F = cellfun(@(x) x(strcmpi(stats{1}(:,1),Factors{fac}),strcmpi(stats{1}(1,:),'F')),stats);F = [F{:}];
                df = cellfun(@(x) x(strcmpi(stats{1}(:,1),Factors{fac}),strcmpi(stats{1}(1,:),'df')),stats);df = [df{:}];
                Clusters{perm,fac} = ClusterExtract(P,SupraTh,F,A);
            end
        else
            for fac = 1:3 % permuting data for each factor and interaction
                % permute factor labels
                clear Fac1p Fac2p;
                switch fac
                    case {1,3}
                        B = arrayfun(@(x) randperm(size(Fac1,2)),1:size(Fac1,1),'uni',false); B = cat(1,B{:});
                        Fac1p = repmat(B,[1 1 size(Fac1,3)]);
                    case {2,3}
                        B = arrayfun(@(x) randperm(size(Fac2,3)),1:size(Fac2,1),'uni',false); B = cat(1,B{:});
                        Fac2p = permute(repmat(B,[1 1 size(Fac2,2)]),[1 3 2]);
                end
                if ~exist('Fac1p','var'), Fac1p = Fac1;end
                if ~exist('Fac2p','var'), Fac2p = Fac2;end    
                % permutation for interaction should be tested later... The
                % best way would be to test interactions on residuals, to
                % avoid the main effects
                % Check it here: http://www.uvm.edu/~dhowell/StatPages/Permutation%20Anova/PermTestsAnova.html
                
                % calculate ANOVA
                for El = 1:size(Data,2)
                    Y = squeeze(Data(:,El,:,:));
                    stats{El} = rm_anova2(Y(:),SubID(:),Fac1p(:),Fac2p(:),Factors(1:2));
                end
                % extract cluster statistics
                P = cellfun(@(x) x(strcmpi(stats{1}(:,1),Factors{fac}),strcmpi(stats{1}(1,:),'P')),stats);P = [P{:}];
                F = cellfun(@(x) x(strcmpi(stats{1}(:,1),Factors{fac}),strcmpi(stats{1}(1,:),'F')),stats);F = [F{:}];
                df = cellfun(@(x) x(strcmpi(stats{1}(:,1),Factors{fac}),strcmpi(stats{1}(1,:),'df')),stats);df = [df{:}];
                C = ClusterExtract(P,SupraTh,F,A);
                Clusters{perm,fac} = C(1); % only the largest cluster
            end
        end
    end
%% make the random distribution and calculate p-vlues for clusters    
    
end