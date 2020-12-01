function [se,zval,pval] = msepd(k1,k2,n1,n2,tails,direction)
% function [se,zval,pval] = msepd(p1,p2,n1,n2,tails,direction)
% 
% computes the standard error of the difference of two proportions (se)
% 
% INPUT
% k1,k2         number of successful outcomes in the two samples; proportions p1 and p2 are therefore k1/n1 and k2/n2
% n1,n2         sample sizes
% direction   
% 
% OUTPUT
% se            standard error of the difference of two proportions
% zval          z-value of the observed difference of proportions
% pval          p-value corresponding to the z-value
% tails         should be 1 or 2 for one- or two-sided p-values
% direction     optional argument, only used when tails==1
%               'p1>p2' and 'p2<p1' allocates the rejection region to the upper range of outcomes
%               'p2>p1' and 'p1<p2' allocates the rejection region to the range of outcomes
% 
% Formulas are from http://stattrek.com/hypothesis-test/difference-in-proportions.aspx on April 4, 2018.
% Calculations crosscheck with http://www.vassarstats.net/propdiff_ind.html on the same day.
% 
% The test statistic (zval) is simple: proportion difference divided by standard error.
% zval = ((p1-p2)-0) / sqrt(p*(1-p)*(1/n1 + 1/n2))
% See http://www.socscistatistics.com/tests/ztest/.
% 
% Maik C. Stuettgen, April 2018
%% input check
if rem(k1,1) ~=0 || rem(k2,1)~=0 || rem(n1,1)~=0 || rem(n2,1)~=0
  error('arguments ''k1'', ''k2'', ''n1'', and ''n2'' must all be integers')
elseif k1>n1
  error('argument ''k1'' must be larger than ''n1''')
elseif k2>n2
  error('argument ''k2'' must be larger than ''n2''')
elseif ~ismember(tails,[1 2])
  error('argument ''tails'' must be either 1 or 2')
elseif exist('direction','var') && all([~strcmp(direction,'p1>p2'),~strcmp(direction,'p1<p2'),~strcmp(direction,'p2>p1'),~strcmp(direction,'p2<p1')])
  error('argument ''direction'' must read ''p1>p2'', ''p2<p1'', ''p1<p2'', or ''p2>p1''')
elseif tails==1 && ~exist('direction','var')
  error('argument ''direction'' not specified')
end
%% the works
% get the pooled sample proportion pp
p1 = k1/n1;
p2 = k2/n2;
pp = (p1*n1 + p2*n2) / (n1+n2);

% compute the standard error of the sampling distribution difference between two proportions
se = sqrt(pp*(1-pp)*(1/n1+1/n2));

% the test statistic is a z-score defined by the following equation:
zval = (p1-p2) / se;

% computation of pval could be shortened, but I like it this way
if tails==1
  if p1>p2 && (strcmp(direction,'p1>p2') || strcmp(direction,'p2<p1')) % results in predicted direction p1>p2
    pval = 1-normcdf(zval);
  elseif p1>p2 && (strcmp(direction,'p1<p2') || strcmp(direction,'p2>p1'))   % results NOT in predicted direction p1>p2
    pval = normcdf(zval);
  elseif p1<p2 && (strcmp(direction,'p1<p2') || strcmp(direction,'p2>p1'))   % results in predicted direction p1<p2
    pval = normcdf(zval);
  elseif p1<p2 && (strcmp(direction,'p1>p2') || strcmp(direction,'p2<p1'))     % results not in predicted direction p1<p2
    pval = 1-normcdf(zval);
  elseif p1==p2
    pval = normcdf(zval);
  end
elseif tails==2
  if p1>p2
    pval = 2*(1-normcdf(zval));
  elseif p1<=p2
    pval = 2*normcdf(zval);
  end
end
end