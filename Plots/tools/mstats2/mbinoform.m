function [f] = mbinoform(n,k,p)
% function [f] = mbinoform(n,k,p)
% 
% binomial formula for computing the exact probability of k occurrences of event A in n trials
% 
% INPUTS:
% n 	number of trials
% k   number of occurrences of event A (can be scalar or vector)
% p   probability of occurrence of event A
% 
% OUTPUT:
% f   probability of observing k occurrences of event A in n trials when the probability of 
%     occurrence of event A equals p
% 
% EXAMPLE:
% 1) Assume have you have a fair coin, and you toss it ten times. What is the probability that
%    heads will come up exactly three times?
%    f = mbinoform(10,3,0.5);
% 2) Same scenario, but this time you are looking for the whole distribution (i.e. 0 to 10
%    occurrences of heads):
%    f = mbinoform(10,0:10,0.5);
%    bar(0:10,f)
% 3) Assume you do not trust that coin. You perform the above experiment, and heads comes up 8 times.
%    Would you consider this statistically significant?
%    You could either use sum(f(9:11)); or 1-binocdf(7,10,0.5)
%    In both cases, the result is 0.0546875 and would thus (just) not be considered
%    "statistically significant" (one-tailed test).
% 
% Maik C. Stuettgen, May 2014
%% input check
if n<k == 1
    error('error - n<k')
end
%% the works
f = nan(numel(k),1);
for i = 1:numel(k)
  f(i) = (factorial(n) / (factorial(k(i))*factorial(n-k(i)))) * p^k(i) * (1-p)^(n-k(i));
end
