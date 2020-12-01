function [sep] = msep(p,n)
% function [sep] = msep(p,n)
% 
% computes the standard error for proportions p for sample size n
% this should not be confused with the standard deviation of the binomial distribution
% the variance of the binomial distribution is n*p*(1-p) (for a series of Bernoulli trials)
% 
% if p and n are vectors of the same size, sep will be computed for p(1) and n(1), p(2) and n(2) etc.
% 
% INPUT
% p     proportion
% n     total number of samples
% 
% Maik C. Stuettgen, July 2014
%% the (little) works
if numel(p)~=numel(n)
  error('Vector inputs do not match in size.')
end
sep = sqrt((p.*(1-p))./n);
end
