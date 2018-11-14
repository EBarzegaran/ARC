function PHStatResults  = RmAnovaPostHoc(Data,StatResults)
% this function performs posthoc analysis on the Data with the same format as
% input to RMAnovaPermute function and StatResults is the output of RMAnovaPermute
% function

% Author: Elham Barzegaran, 11/14/2018
%%
for f = 1:numel(StatResults) % for each main effect and interaction
    Stats = StatResults{f};
    for cl = 1:numel(Stats.Clusters)
        if Stats.Clusters(cl).Pvalue<0.05
            SE = Stats.Clusters.Nodes;
            for El = 1:numel(SE)
                % prepare Data and fit the model usinf fit rm
                Y = squeeze(Data(:,SE(El),:,:));
                varNames = arrayfun(@(x) ['y' num2str(x)],1:(size(Y,2)*size(Y,3)),'uni',false);
                t = array2table(reshape(Y,size(Y,1),(size(Y,2)*size(Y,3))),'VariableNames',varNames);
                % prepare factors
                factorNames = {'fac1','fac2'};
                L1s = repmat(1:size(Y,2),1,size(Y,3));
                L1s = arrayfun(@(x) ['L1_' num2str(x)],L1s,'uni',false);
                L2s = reshape(repmat(1:size(Y,3),size(Y,2),1),1,(size(Y,2)*size(Y,3)));
                L2s = arrayfun(@(x) ['L2_' num2str(x)],L2s,'uni',false);
                within = table(L1s',L2s','VariableNames',factorNames);
                rm = fitrm(t,['y1-y' num2str(size(Y,2)*size(Y,3)) '~1'],'WithinDesign',within);
                % this have the same result as rmanova2
                %[ranovatbl] = ranova(rm, 'WithinModel','fac1+fac2+fac1*fac2');
                
                switch Stats.Factor
                    case {'fac1','fac2'}
                        T{El} = multcompare(rm,Stats.Factor);
                    case 'fac1 x fac2'
                        T{El,1} = multcompare(rm,'fac2','By','fac1');
                        T{El,2} = multcompare(rm,'fac1','By','fac2');
                    otherwise
                        error('Unrecognize factor for posthoc analysis');
                end
            end
            StatResults{f}.Clusters(cl).Posthoc = T;
            clear T;
        end
    end
end
end