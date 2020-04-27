function StatResults  = RmAnovaPostHoc(Data,StatResults)
% this function performs posthoc analysis on the Data with the same format as
% input to RMAnovaPermute function and StatResults is the output of RMAnovaPermute
% function

% Author: Elham Barzegaran, 11/14/2018
%%
for f = 1:numel(StatResults) % for each main effect and interaction
    Stats = StatResults{f};
    for cl = 1:numel(Stats.Clusters) % for each cluster
        Num=1;
        if Stats.Clusters(cl).Pvalue<=0.05
            SE = Stats.Clusters(cl).Nodes;
            Y = squeeze(Data(:,SE,:,:));% should we correct for all electrodes or just the ones that survive?
            switch Stats.Factor 
                 case {'fac1','fac2'} % main effect post-hoc
                     FN = str2double(Stats.Factor(4));% Factor Number
                     if size(Y,FN+2)>1
                        Kfl = nchoosek(1:size(Y,FN+2),2);
                        for l = 1:size(Kfl,1) % do every pair of comparision
                            inds = arrayfun(@(x) 1:size(Y,x),1:ndims(Y),'uni',false);
                            inds{FN+2} = Kfl(l,:);                           
                            %stats = RmAnovaPermute(Y(inds{:}),A(SE,SE),PermNum,SupraTh,StatType,FN);
                            results = meanrmanova(Y(inds{:}),Stats.Factor);
                            results.levels = Kfl(l,:); 
                            results.factor = FN;
                            Stats.Clusters(cl).PostHoc(l) = results;
                        end
                     end
                  case 'fac1 x fac2' % interaction post-hocs
                      for FN1 = 1:2 % over factors
                          if size(Y,FN1+2)>1 
                              FN2 = setdiff(1:2,FN1);
                              for f1l = 1:size(Y,FN1+2)                                 
                                  Kf2l = nchoosek(1:size(Y,FN2+2),2);
                                  for l2 = 1:size(Kf2l,1)% do every pair of comparision
                                      inds = arrayfun(@(x) 1:size(Y,x),1:ndims(Y),'uni',false);
                                      inds{FN1+2} = f1l;
                                      inds{FN2+2} = Kf2l(l2,:);
                                      %PHStatResults = RmAnovaPermute(Y(inds{:}),A(SE,SE),PermNum,SupraTh,StatType,FN2);
                                      results = meanrmanova(Y(inds{:}),['fac' num2str(FN2)]);
                                      results.factor1 = FN1;
                                      results.levels1 = f1l;
                                      results.factor2 = FN2;
                                      results.levels2 = Kf2l(l2,:);
                                      Stats.Clusters(cl).PostHoc(Num)= results;
                                      Num = Num+1;
                                  end
                              end
                          end
                      end
                      
                otherwise
                    error('Unrecognize factor for posthoc analysis');
            end

        end
    end
    StatResults{f}=Stats;
end

end

function results = meanrmanova(y,FactorName)
for el = 1:size(y,2)
    Y = squeeze(y(:,el,:,:));Y = reshape(Y,size(y,1),size(y,3),size(y,4));
    % mean calculations
    factornum = str2double(FactorName(4));
    Y2 = permute(Y,[ setdiff(1:ndims(Y),factornum+1) factornum+1]);% first dimension is subjects\
    Y2 = reshape(Y2,size(Y2,2)*size(Y2,1),size(Y2,3));
    results.mean = mean(Y2); % mean of data at different levels
    results.SEM = std(Y2)/sqrt(size(Y2,1));% standard error of mean
    
    % prepare for anova
    varNames = arrayfun(@(x) ['y' num2str(x)],1:(size(Y,2)*size(Y,3)),'uni',false);
    t = array2table(reshape(Y,size(Y,1),(size(Y,2)*size(Y,3))),'VariableNames',varNames);
    factorNames = {'fac1','fac2'};
    L1s = repmat(1:size(Y,2),1,size(Y,3));
    L1s = arrayfun(@(x) ['L1_' num2str(x)],L1s,'uni',false);
    L2s = reshape(repmat(1:size(Y,3),size(Y,2),1),1,(size(Y,2)*size(Y,3)));
    L2s = arrayfun(@(x) ['L2_' num2str(x)],L2s,'uni',false);
    within = table(L1s',L2s','VariableNames',factorNames);
    rm = fitrm(t,['y1-y' num2str(size(Y,2)*size(Y,3)) '~1'],'WithinDesign',within);
    [ranovatbl] = ranova(rm, 'WithinModel',FactorName);
    results.F(el) = ranovatbl.F(3);
    results.P(el) = ranovatbl.pValue(3);
    results.df(el,:) = ranovatbl.DF(2:3);
end
%              
end