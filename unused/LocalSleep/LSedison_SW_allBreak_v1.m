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

data_path='/Users/tand0009/Data/LS_Edison/EDF_fixed';
save_path='/Users/tand0009/Data/LS_Edison/LocalSleep';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*.edf']);

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).


%%
all_SW=[];
all_Subds=[];
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
    load([Folder_Name filesep 'T_' SubdID '.mat'])
    %     scores=find(T.StadeL~='?');
    fprintf('... processing %s\n',File_Name);
    
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([Folder_Name filesep File_Name]);
    events=ft_read_event([Folder_Name filesep File_Name]);
    %         dat=ft_read_data([Folder_Name filesep File_Name]);
    
    %%% consider only the pause
    %     hypnogram = import_hypnogram([Folder_Name filesep File_Name(1:4) '_hypnogram.TXT']);
    
    load([save_path filesep File_Name(1:end-4) '_allSW'],'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=9; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=100;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    if paramSW.AmpCriterionIdx==9
    all_Waves(:,9)=abs(all_Waves(:,9));
    end
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:3
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    
    start_pause=(T.Epoch(find(T.Start(:,1)))-1)*30;
    end_pause=(T.Epoch(find(T.End(:,1))))*30;
    
    
    all_SW=[all_SW ; [sum(slow_Waves(:,3)==3)/(end_pause-start_pause)*60 nanmean(slow_Waves(slow_Waves(:,3)==3,4)) nanmean(slow_Waves(slow_Waves(:,3)==3,12)) nanmean(slow_Waves(slow_Waves(:,3)==3,13))]];
    all_Subds=[all_Subds ; {SubdID}];
    
end


%%
data_clean=readtable('../Edison_Tables/Clean_Data_SS_Pow.txt');
data_clean.LSdens=nan(size(data_clean,1),1);
data_clean.LSamp=nan(size(data_clean,1),1);
data_clean.LSdslp=nan(size(data_clean,1),1);
data_clean.LSuslp=nan(size(data_clean,1),1);

for nS=1:length(all_Subds)
    this_line=match_str(data_clean.ID,all_Subds(nS));
    
    data_clean.LSdens(this_line)=all_SW(nS,1);
    data_clean.LSamp(this_line)=all_SW(nS,2);
    data_clean.LSdslp(this_line)=all_SW(nS,3);
    data_clean.LSuslp(this_line)=all_SW(nS,4);
    
    
end

writetable(data_clean,'../Edison_Tables/Clean_Data_SS_Pow_LS.txt');