# Edison_scripts

The following are the primary scripts that were used to analyze the data and generate the figures in the study "Sleep onset is a creative sweet spot"
All raw data (EEG and behavioral) are available at the following link:

There are several main steps: 

1) You can generate a T-file for each subject with the folder "Pack-Generate_EEG_Tables" and the script "Generate_EEG_Tables".
This is a 1-by-1s table that contains all of the information of the nap (start and end timings, sleep scoring, bottle drops...). 
You can skip this step by directly downloading the T files straight from OSF.

2) If you run both scripts "Analyses_Edison_EEG" and "Anayses_Edison_Behaviour", you will generate tables with all of the necessary information needed to 
generate behavioral figures as well as demographic data. You can also skip this step by directly downloading this table ("All_Data.mat") on OSF.

All variable names should be self-explanatory, but here is a quick rundown of some of the more essential variables in this table:

Data.InsightPre and Data.InsightPost indicate whether or not subjects had insight (pre and post).
Data.AhaMoment corresponds to the Eureka moment (detected by the matlab function "findchangepts").
Data.RT_beforeAha and Data.RT_fromAha correspond to the mean RT obtained before and after the Eureka moment, respectively.
Data.NbEpochsN1Edison and Data.NbEpochsN2Edison give you the number of epochs obtained in each subtage.
Subjects with one or more epochs of N1 are assigned to the N1 group, subjects with one or one epoch of N2 (with or without N1 epochs) are assigned to the N2 group,
and the remaining subjects are assigned to the Wake group.
Data.Bottle indicates whether or not the subjects dropped the bottle.
Data.Hypna, Data.HypnaTask and Data.HypnaBottle correspond to whether they reported an hypnagogia, or a task-related hypnagogia, respectively.
Data.Metasleep indicates whether the subjects subjectively perceived themselves to be asleep when the bottle dropped.


3) EEG Analyses




