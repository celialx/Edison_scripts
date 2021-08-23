# Edison_scripts

This README presents the primary scripts that were used to analyze the data and generate the figures in the study "Sleep onset is a creative sweet spot"
These scripts are located in the GitHub repository "Edison_scripts" () in the subfolder Core_Scripts_Article"
All raw data (EEG and behavioral) are available at the following link:
https://osf.io/dkbw7

There are several main steps: 

1) You can generate a T-file for each subject using the script "Generate_EEG_Tables" and the subfolder with the folder "Pack-Generate_EEG_Tables"

This script creates a table that contains all of the information of the nap of each participant (start and end timings, sleep scoring, bottle drops...). 
You can skip this step by directly downloading the T files from OSF (component "T matrix (scoring)").

2) If you run both scripts "Analyses_Edison_EEG_Clean" and "Anayses_Edison_Behaviour_Clean", you will generate tables with all of the necessary information needed to 
generate behavioral figures as well as demographic data. You can also skip this step by directly downloading this table ("All_Data_Behaviour.mat") on OSF (component "Tables").

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


3) EEG Analyses are composed of two main analyses: an analysis on the continuous data of each nap (Figure 3) and of the data locked on bottle drops (Figure 4).

Elements of Figure 3 are generated using the function "LSedison_TF_byEpoch_vFigure_submitted.m" for Matlab.
Elements of Figure 4 are generated using the function "LSedison_TF_aroundDrop_vFigure2_submitted" for Matlab. 

Dependencies include:
- the Fiedltrip toolbox (https://github.com/fieldtrip/fieldtrip),
- the FOOOF and FOOOF-mat toolboxes (https://github.com/fooof-tools/fooof and https://github.com/fooof-tools/fooof_mat),
- the LSCPtools toolbox (https://github.com/andrillon/LSCPtools) and,
- the export_fig (https://github.com/altmany/export_fig) toolbox.

These scripts run from the EDF files (available on OSF in the component "Raw EEG (EDF)"), the T matrices (see point 1 and component "T matrix (scoring)") and the behavioural table (component "Tables"). Path will need to be edited to match your local machine.  

Statistics and models are performed by the script "LSedison_LMEs_Pow_v2_submitted.Rmd" for R. Dependencies are listed at the beginning of the script. The script uses the table Clean_Data_SS_Pow_submitted.txt (created by LSedison_TF_byEpoch_vFigure_submitted.m but available on OSF, component "Tables"). 





