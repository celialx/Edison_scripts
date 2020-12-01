function [eta2,omega2] = meta2(values,groups)
% function [eta2,omega2] = meta2(counts,groups)
%	computes eta squared OR omega squared for n groups
% just to be quicker than mes1way
% 
% values      vector which contains values of the dependent variable (real numbers)
% groups      vector which contains group assignments (integers)
% 
% eta2 = SSbetween/SStotal = SSB/SST
% sst = sum((xij-mean(x))^2)
% ssb = sum(ni*(mean(xi)-mean(x))^2)
% 
% omega2 = (df_effect*(MS_effect - MS_W)) / (SStotal + MS_W)
% formula for omega2 taken from Kline, Beyond significance testing; result identical to mes1way
% 
% Maik C. Stuettgen
%% the works
sst      = sum((values-nanmean(values)).^2);
groupIDs = unique(groups);
ssgroups = nan(numel(groupIDs),1);
for i = 1:numel(ssgroups)
  ssgroups(i) = sum(groups==groupIDs(i)) * (nanmean(values(groups==groupIDs(i))) - nanmean(values)).^2;
end
ssb  = sum(ssgroups);
eta2 = ssb/sst;
[~,tbl] = anova1(values,groups,'off');
omega2 = (tbl{2,3}*(tbl{2,4}-tbl{3,4})) / (tbl{4,2}+tbl{3,4}); % formula from Rex Kline

% the formulas below all yield the same value:
% formula 28.10 from Diehl and Arbinger, p. 651
% omega2b = (tbl{2,2} - (numel(ssgroups)-1)*tbl{3,4}) / (tbl{4,2}+tbl{3,4});
% formula 28.11 from Diehl and Arbinger, p. 651
% omega2c = ((numel(ssgroups)-1)*(tbl{2,5}-1)) / ((numel(ssgroups)-1)*(tbl{2,5}-1)+numel(values));
% according to Diehl & Arbinger (p. 652, formula 28.17), omega2 can be readily computed from eta2
% omega2d = (eta2*(numel(values)-1) - (numel(ssgroups)-1)) / (1-eta2+(numel(values)-numel(ssgroups)));