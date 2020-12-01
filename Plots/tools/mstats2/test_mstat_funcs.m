%% test mstat functions
%% binocof
x = nchoosek(10,4);
y = binocof(10,4);
%% binoform
% http://www.mathwords.com/b/binomial_probability_formula.htm
% You are taking a 10 question multiple choice test. If each question has four choices 
% and you guess on each question, what is the probability of getting exactly 7 questions correct?
% p = 0.0031
p = binoform(10,7,0.25);
%% cells2vectors
load mstats_exampledata
for i = 1:3
  disp(x{i})
end
y = cell2vectors(x);
anova1(y(:,1),y(:,2))
%% fishertest
% http://en.wikipedia.org/wiki/Fisher%27s_exact_test
% A sample of teenagers might be divided into male and female on the one hand, and those 
% that are and are not currently dieting on the other. We hypothesize that the proportion 
% of dieting individuals is higher among the women than among the men, and we want to test 
% whether any difference of proportions that we observe is significant.
%           MEN WOM
%          +---+---+   
% DIET     | a | b |
%          +---+---+
% NO DIET  | c | d |
%          +---+---+  
% a=1;b=9;c=11;d=3;
% p = fishertest([1,9;11,3]);
% p = 0.0013461
%% rci
% example from Diehl and Arbinger, Einführung in die Inferenzstatistik (2nd edition), pp. 376
[cilohi,p] = rci(0.41,648,0.01,0);
%% rddiffci
% example from Bortz, Statistik (4th edition), pp. 204
% (note that Bortz provides z values only)
r13   = 0.41;
r23   = 0.52;
r12   = 0.48;
n     = 100;
alpha = 0.05;
[rddiff,cilohi,p] = rddiffci(r13,r23,r12,n,alpha);

% example from Diehl and Arbinger, Einführung in die Inferenzstatistik (2nd edition), pp. 382
% (note that the authors provide t values only)
[rddiff,cilohi,p] = rddiffci(-0.336,-0.126,0.413,342,alpha);

% example from Zou paper
% (author provides values for another method, but these are quite close)
[rddiff,cilohi,p] = rddiffci(.396,.179,0.088,66,0.05);
%% ridiffci
% example from Bortz, Statistik (4th edition), pp. 203
% author provides critical z value only
r1 = .38;
r2 = .65;
n1 = 60;
n2 = 40;
alpha = 0.05;
[rdiff,cilohi,p] = ridiffci(r1,r2,n1,n2,alpha);