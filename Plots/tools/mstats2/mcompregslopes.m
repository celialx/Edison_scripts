function [p,stats] = mcompregslopes(x1,y1,x2,y2)
% function [p,stats] = mcompregslopes(x1,y1,x2,y2)
% Function compares the slopes of two regression lines.
%
% INPUT
% x1,y1,x2,y2   column vectors in which numel(x1)=numel(y1) and numel(x2)=numel(y2)
%               regression will be computed with y as criterion and x as predictor variable
%
% OUTPUT
% p       p-value of the comparison
% stats   struct containing fields for the corresponding t-value and degrees of freedom (.t and .df, where df=numel(x1)+numel(x2)-4
%
% validated with data taken from:
% http://www.real-statistics.com/regression/hypothesis-testing-significance-regression-line-slope/comparing-slopes-two-independent-samples/
% 
% Maik C. Stuettgen, March 2016
%% input check
if numel(x1)~=numel(y1)
  error('x1 and y1 have unequal numbers of elements')
elseif numel(x2)~=numel(y2)
  error('x2 and y2 have unequal numbers of elements')
end
if size(x1,1)<size(x1,2),x1=x1';end
if size(y1,1)<size(y1,2),y1=y1';end
if size(x2,1)<size(x2,2),x2=x2';end
if size(y2,1)<size(y2,2),y2=y2';end
if ~isvector(x1) || ~isvector(y1) || ~isvector(x2) || ~isvector(y2)
  error('At least one of the input arrays is not a vector.')
end
%% compute regression coefficients
a1b1 = polyfit(x1,y1,1);
a2b2 = polyfit(x2,y2,1);
disp(['B_x1y1 = ',num2str(a1b1(1),'%2.2f'),', intercept = ',num2str(a1b1(2),'%2.2f')])
disp(['B_x2y2 = ',num2str(a2b2(1),'%2.2f'),', intercept = ',num2str(a2b2(2),'%2.2f')])
%% compare regression slopes
stats.t  = (a1b1(1)-a2b2(1)) / sqrt(sebyx(x1,y1)^2 + sebyx(x2,y2)^2);  % when testing the differences, variances add up
stats.df = numel(x1)+numel(x2)-4;
p = 2*(1-tcdf(abs(stats.t),stats.df));                   % multiply by 2 because we employ a two-tailed test
disp(['t(',num2str(stats.df),') = ',num2str(stats.t,'%2.2f')])
disp(['p = ',num2str(p,'%2.5f')])
end


%% HELPER FUNCTION -------------------------------------------------------------------------------------------
function se = sebyx(x,y)
% compute the standard error of a regression coefficient
if size(x,1)<size(x,2) || size(y,1)<size(y,2)
  % x and y need to be column vectors
  error('provide column vectors')
end
se = (std(y) / std(x)) * sqrt((1-corr(x,y)^2)/(numel(x)-2));
end