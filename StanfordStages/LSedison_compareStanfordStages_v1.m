%%
clear all
close all

path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
path_localsleep='/Users/tand0009/WorkGit/projects/inprogress/wanderIM/localsleep';
addpath(path_fieldtrip);
addpath(path_localsleep);
ft_defaults;

path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

data_path='/Users/tand0009/Data/LS_Edison/';
save_path='/Users/tand0009/Data/LS_Edison/LocalSleep';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*' filesep '*.edf']);

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'int�gralit� de l'exp�rience donc la pause ne commence normalement pas au d�but de l'enregistrement (sauf si on a oubli� de lancer l'enregistrement au tout d�but...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes repr�sente le d�but et la fin de la pause. Tu souhaites donc regarder les trac�s EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces m�mes infos de mani�re peut-�tre plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).


%%
doERP=1;
for nF=1:length(files)
    File_Name=files(nF).name;
    Folder_Name=files(nF).folder;
    SubdID=File_Name;
    sep=findstr(File_Name,'_');
    SubdID=SubdID(1:sep(1)-1);
    
    if exist([Folder_Name filesep 'T_' SubdID '.mat'])==0
        warning('matrix file missing');
        continue;
    end
    if exist([Folder_Name filesep '..' filesep 'StanfordStages' filesep File_Name(1:end-4) '.hypnogram.txt'])==0
        warning('Stanford-stages file missing');
        continue;
    end
    load([Folder_Name filesep 'T_' SubdID '.mat'])
    hypnodensity = import_stage_hypnodensity([Folder_Name filesep '..' filesep 'StanfordStages' filesep File_Name(1:end-4) '.hypnodensity.txt']);
    hypnogram = import_stage_hypnogram([Folder_Name filesep '..' filesep 'StanfordStages' filesep File_Name(1:end-4) '.hypnogram.txt']);
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([Folder_Name filesep File_Name]);
    
    figure;
    ASSM_Score=T.Stade;
    Beg_Task=T.TimeID(find(T.Start(:,1)));
    End_Task=T.TimeID(find(T.End(:,1)));
    ASSM_Score=ASSM_Score(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    ASSM_Score_Time=T.TimeID(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    
    plot(ASSM_Score_Time,ASSM_Score);
    hold on;
    hypnogram_perSec=reshape(repmat(hypnogram,1,15)',1,numel(repmat(hypnogram,1,15)));
    time_hypno_perSec=1:length(hypnogram_perSec);
    Stanford_Score=hypnogram_perSec(time_hypno_perSec>=Beg_Task & time_hypno_perSec<=End_Task);
    Stanford_Score_Time=time_hypno_perSec(time_hypno_perSec>=Beg_Task & time_hypno_perSec<=End_Task);
    
    
    plot(Stanford_Score_Time,Stanford_Score,'LineStyle','--','LineWidth',2);
    xlabel('Time (s)')
    ylabel('Stage')
    ylim([-0.5 5.5])
    set(gca,'YTick',[0 1 2 3 5],'YTickLabel',{'W','N1','N2','N3','REM'});
    format_fig;
    title(SubdID)
    legend('manual','stanford')
end