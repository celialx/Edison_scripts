%%
clear all
close all

%% Test Importance of Sleep and co-variates
load ../Edison_Tables/All_Data.mat

data_clean=Data(Data.InsightPre==0,:);
data_clean(:,[20 41 62 81 ])=[];
% writetable(data_clean,'Edison_Tables/Clean_Data.txt');

%%
mdl0=fitglm(data_clean,'InsightPost~1+SleepEdison+Age+Gender+Epworth+DreamRecall+EduLevel+Laterality+UsedToEnigmas')

%%
mdl1=fitglm(data_clean,'InsightPost~1+SleepEdison+Bottle+Age+Gender+Epworth+DreamRecall+EduLevel+Laterality+UsedToEnigmas','Distribution','binomial')
