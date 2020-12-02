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
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).


%%
doERP=1;
all_TF_drops=[];
StageDrop=[];
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
    if ~isempty(find(T.Ball(:,1)))
        fprintf('... processing %s\n',File_Name);
        
        
        % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
        hdr=ft_read_header([Folder_Name filesep File_Name]);
        events=ft_read_event([Folder_Name filesep File_Name]);
        %         dat=ft_read_data([Folder_Name filesep File_Name]);
        
        %%% consider only the pause
        %     hypnogram = import_hypnogram([Folder_Name filesep File_Name(1:4) '_hypnogram.TXT']);
        
        cfg=[];
        cfg.trialfun             = 'LSedison_trialfun_drops';
        cfg.dataset             = [Folder_Name filesep  File_Name];
        cfg.trialdef.prestim    = 2*60;
        cfg.trialdef.poststim   = .5*60;
        cfg.demean          = 'yes';
        cfg.T               = T;
        cfg = ft_definetrial(cfg);
        cfg.channel={'Fp1-A2','C3-A2','O1-A2'};
        data                   = ft_preprocessing(cfg); % read raw data
        
        %          trl = LSedison_trialfun_drops(cfg)
        %         start_drop=(T.TimeID(find(T.Ball(:,1)))-1)*hdr.Fs;
        
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
        cfg.toi          = [-2*60:2:0.5*60];                         % time
        cfg.keeptrials  = 'yes';
        TFRhann = ft_freqanalysis(cfg, data);
        
        temp=(log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(:,:,:,:),3),[1 1 size(TFRhann.powspctrm,3) 1])));
        
        all_TF_drops=cat(1,all_TF_drops,temp);
        
        stages_drop=T.Stade(find(T.Ball(:,1)));
        for k=1:length(stages_drop)
            StageDrop=[StageDrop stages_drop(k)];
        end
    end
    %     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
end

%%
freqs=TFRhann.freq;
times=TFRhann.time;
temp_toplot=squeeze(mean(all_TF_drops(:,2,:,:),1));
h=simpleTFplot(temp_toplot,freqs,times,0,1);
caxis([-4 3])
colorbar;

xlabel('Time from Drop (s)')
ylabel('Frequency (Hz)')
title(sprintf('N=%g drops',size(all_TF_drops,1)))
format_fig;
 %%
figure;
subplot(1,2,1);
temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,2,:,:),1));
 h=simpleTFplot(temp_toplot,freqs,times,0,0);
 caxis([-4 4])
 
 subplot(1,2,2);
temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,2,:,:),1));
 h=simpleTFplot(temp_toplot,freqs,times,0,0);
 caxis([-4 4])
 
 %%
 FOI=[1 4];
 COI=2;
 figure;
 temp_toplot=squeeze(mean(all_TF_drops(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
 simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,0,1,1);
 xlim([-60 20])
 format_fig;
 xlabel('Time from Drop (s)')
 ylabel('Power (dB)')
 title(sprintf('N=%g drops',size(all_TF_drops,1)))
%  hold on;
%  temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,0,1,1);
%  temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'b',0,'-',0.5,1,0,1,1);
%  