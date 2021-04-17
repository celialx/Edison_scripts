%%
clear all
close all

path_localsleep='D:\MATLAB\Toolbox\wanderIM\localsleep\';
addpath(path_localsleep);
ft_defaults;

path_LSCPtools='D:\MATLAB\Toolbox\LSCPtools\';
addpath(genpath(path_LSCPtools));

 data_path='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\EEG\';
save_path='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\LocalSleep\'
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path '*.edf']);

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).

%% FOOOF
fooof_path='/Users/tand0009/WorkGit/projects/ext/fooof_mat/';
f_range = [2, 30];
settings = struct();  % Use defaults
addpath(genpath(fooof_path));

%%
doERP=1;
all_fooof_drops=cell(3,2);
all_spec_drops=cell(3,1);
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
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:0.1:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*10;   % length of time window = 0.5 sec
        cfg.toi          = [-2*60:2:0.5*60];                         % time
        cfg.keeptrials  = 'yes';
        TFRhann = ft_freqanalysis(cfg, data);
        
        freqs=TFRhann.freq;
        times=TFRhann.time;
        power=TFRhann.powspctrm;
        for k=1:size(power,1)
            for nEl=1:size(power,2)
                temp_fooof=[];
                for i=1:length(times)
                    temp_power=squeeze(power(k,nEl,:,i));
                    if sum(isnan(temp_power))>0
                        temp_fooof(i,:)=nan(1,2);
                        continue;
                    end
                    try
                        fooof_results = fooof(freqs, temp_power', f_range, settings,1);
                        temp_fooof(i,:)=fooof_results.background_params;
                    catch
                        temp_fooof(i,:)=nan(1,2);
                    end
                end
                all_fooof_drops{nEl,1}=[ all_fooof_drops{nEl,1} ; temp_fooof(:,1)'];
                all_fooof_drops{nEl,2}=[ all_fooof_drops{nEl,2} ; temp_fooof(:,2)'];
                all_spec_drops{nEl}=cat(3,all_spec_drops{nEl},squeeze(nanmean(power(:,nEl,:,:),1)));
            end
        end
        stages_drop=T.Stade(find(T.Ball(:,1)));
        for k=1:length(stages_drop)
            StageDrop=[StageDrop stages_drop(k)];
        end
    end
    
end

%%
figure;
nEl=2;
% for nEl=1:3
%    subplot(1,3,nEl);
format_fig;
temp_plot=all_fooof_drops{nEl,2};
simpleTplot(times,temp_plot,0,'k',0,'-',0.5,1,0,1,1);
xlim([-60 20])
% end


%%
figure;
hold on;
cmap=cbrewer('seq','OrRd',12);
for k=1:12
    plot(freqs,squeeze(nanmean(nanmean(all_spec_drops{nEl}(:,((k-1)*5+1):(5*k),:),2),3)),'Color',cmap(k,:));
end