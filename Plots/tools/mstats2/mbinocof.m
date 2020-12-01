function [bc] = mbinocof(n,k)
% function [bc] = mbinocof(n,k)
% 
% computes the binomial coefficient for values of n up to 170
%
% In combinatorics, (n k)' is interpreted as the number of k-element subsets (the k-combinations) of an n-element
% set, that is the number of ways that k things can be 'chosen' from a set of n things. Hence, (n k)' is often read as
% "n choose k" and is called the choose function of n and k.
%
% Sadly, this code became obsolete when The Mathworks introduced nchoosek.m.
% 
% INPUTS
% n 	number of trials
% k   number of occurrences of event A
%
% OUTPUT
% bc  binomical coefficient
%
% Code validated with Pascal's triangle:
% n           ks from 0 to n
% 0                1
% 1             1     1
% 2           1    2     1
% 3         1    3    3    1
% 4       1   4    6    4    1
% 5     1   5   10   10   5    1
%
% Further validation tests were run in comparison to nchoosek.
% 
% EXAMPLE
% The binomial coefficient is the number of ways of picking k unordered outcomes from n possibilities.
% This can be nicely illustrated with pizza toppings. A local pizza service offers 31 different toppings.
% Suppose you want to have a pizza with six (different!) toppings.
% Type mbinocof(31,6) (or nchoosek(31,6)). The choice will be hard:
% There are n!/(k!(n-k)!) = 31!/(6!(31-6)!) = 736,281 different pizzas to choose from.
% 
% Maik C. Stuettgen, May 2014
%% input check
if n<k
  error('error - n<k!')
elseif mod(n,1)~=0 || mod(k,1)~=0
  error('error - input values are not integers')
end
%% the works
if ~isa(n,'double'),n=double(n);end
if ~isa(k,'double'),k=double(k);end
bc = factorial(n) / (factorial(k)*factorial(n-k));