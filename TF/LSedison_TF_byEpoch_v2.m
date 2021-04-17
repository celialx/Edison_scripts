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

%% FOOOF
fooof_path='/Users/tand0009/WorkGit/projects/ext/fooof_mat/';
f_range = [2, 30];
settings = struct();  % Use defaults
addpath(genpath(fooof_path));


%%
nc=0;
StageDrop=[];
all_Pow=[];
all_Subds=[];
all_ERP=[];
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
    
    cfg=[];
    cfg.trialfun             = 'LSedison_trialfun_break_contepochs';
    cfg.dataset             = [Folder_Name filesep  File_Name];
    cfg.demean          = 'yes';
    cfg.T               = T;
    cfg.lengthepochs = 30;
    cfg = ft_definetrial(cfg);
    cfg.channel={'Fp1-A2','C3-A2','O1-A2'};
    data                   = ft_preprocessing(cfg); % read raw data
    
    maxAbs=[];
    for k=1:length(data.trial)
        maxAbs(k,:)=max(abs(data.trial{k}),[],2);
    end
    %     keep=ones(1,length(data.trial));
    %     for k=1:length(data.trial)
    %         %         all_ERP=[all_ERP ; data.trial{k}(2,:)];
    %
    %         if max(abs(data.trial{k}(2,data.time{k}<0)))>250
    %             keep(k)=0;
    %         end
    %     end
    %     if isempty(find(keep))
    %         continue;
    %     end
    nc=nc+1;
    %          trl = LSedison_trialfun_drops(cfg)
    %         start_drop=(T.TimeID(find(T.Ball(:,1)))-1)*hdr.Fs;
    
    threshArt=100;
    fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(max(maxAbs,[],2)>threshArt),100*mean(max(maxAbs,[],2)>threshArt),threshArt)
    
    cfg              = [];
    cfg.trials       = find(max(maxAbs,[],2)<=threshArt);
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
    cfg.toi          = 0:0.2:30;                         % time
    cfg.keeptrials  = 'yes';
    TFRhann = ft_freqanalysis(cfg, data);
    
    temp=(log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(:,:,:,:),3),[1 1 size(TFRhann.powspctrm,3) 1])));
    
    AASM_Score=T.Stade;
    Beg_Task=(T.Epoch(find(T.Start(:,1))));
    End_Task=T.Epoch(find(T.End(:,1)));
    AASM_Score=AASM_Score(T.Epoch>=Beg_Task & T.Epoch<=End_Task);
    AASM_Score_Epochs=(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score_Epochs2=unique(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score2=nan(1,length(AASM_Score_Epochs2));
    for j=1:length(AASM_Score_Epochs2)
        if sum(isnan(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)))))==length(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j))))
            AASM_Score2(j)=nan;
        else
            AASM_Score2(j)=unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)));
        end
    end
    AASM_Score2=AASM_Score2(1:size(temp,4));
    freqs=TFRhann.freq;
    %         for j=1:length(AASM_Score2)
    %             temp_power=squeeze(TFRhann.powspctrm(1,2,:,j));
    %             freqs=TFRhann.freq;
    %             fooof_results = fooof(freqs, temp_power', f_range, settings,1);
    %
    %         end
    Stages=0:2;
    temp_power=squeeze(TFRhann.powspctrm(1,2,:,:));
    fooof_results_perStage=[];
    for nSta=1:3
        Pow_perStage(nc,nSta,:)=squeeze(nanmedian(temp_power(:,AASM_Score2==Stages(nSta)),2));
        if sum(isnan(squeeze(Pow_perStage(nc,nSta,:))'))==length(squeeze(Pow_perStage(nc,nSta,:))')
            fooof_results_perStage{nSta}.background_params=nan(1,2);
        else
            fooof_results_perStage{nSta} = fooof(freqs,squeeze(Pow_perStage(nc,nSta,:))', f_range, settings,1);
        end
    end
    Pow_all(nc,:)=log(squeeze(nanmean(temp_power,2)));
    fooof_results = fooof(freqs,squeeze(nanmedian(temp_power,2))', f_range, settings,1);
    
    all_Pow=[all_Pow ; [nF median(Pow_all(nc,freqs>1 & freqs<4)) median(Pow_all(nc,freqs>4 & freqs<8)) median(Pow_all(nc,freqs>8 & freqs<11)) median(Pow_all(nc,freqs>1 & freqs<40)) fooof_results.background_params length(AASM_Score2) sum(AASM_Score2==0) sum(AASM_Score2==1) sum(AASM_Score2==2) ...
        median(log(Pow_perStage(nc,1,freqs>1 & freqs<4))) median(log(Pow_perStage(nc,1,freqs>4 & freqs<8))) median(log(Pow_perStage(nc,1,freqs>8 & freqs<11))) median(log(Pow_perStage(nc,1,freqs>1 & freqs<40))) fooof_results_perStage{1}.background_params...
        median(log(Pow_perStage(nc,2,freqs>1 & freqs<4))) median(log(Pow_perStage(nc,2,freqs>4 & freqs<8))) median(log(Pow_perStage(nc,2,freqs>8 & freqs<11))) median(log(Pow_perStage(nc,2,freqs>1 & freqs<40))) fooof_results_perStage{2}.background_params...
        median(log(Pow_perStage(nc,3,freqs>1 & freqs<4))) median(log(Pow_perStage(nc,3,freqs>4 & freqs<8))) median(log(Pow_perStage(nc,3,freqs>8 & freqs<11))) median(log(Pow_perStage(nc,3,freqs>1 & freqs<40))) fooof_results_perStage{3}.background_params ]];
    all_Subds=[all_Subds ; {SubdID}];
    %     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
end

%%
freqs=TFRhann.freq;
times=TFRhann.time;
figure; hold on;
for nSta=1:3
    temp_toplot=squeeze(Pow_perStage(:,nSta,:));
    temp_toplot=log10(temp_toplot./repmat(nanmean(temp_toplot,2),1,length(freqs)));
    plot(freqs,nanmean(temp_toplot));
end
xlabel('Frequency (Hz)')
format_fig;


%%
data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');
data_clean.PowDelta=nan(size(data_clean,1),1);
data_clean.PowTheta=nan(size(data_clean,1),1);
data_clean.PowAlpha=nan(size(data_clean,1),1);

data_clean.PowDeltaN=nan(size(data_clean,1),1);
data_clean.PowThetaN=nan(size(data_clean,1),1);
data_clean.PowAlphaN=nan(size(data_clean,1),1);

data_clean.PowDeltaN_W=nan(size(data_clean,1),1);
data_clean.PowThetaN_W=nan(size(data_clean,1),1);
data_clean.PowAlphaN_W=nan(size(data_clean,1),1);

data_clean.PowSlope=nan(size(data_clean,1),1);
data_clean.PowBG=nan(size(data_clean,1),1);

for nS=1:length(all_Subds)
    this_line=match_str(data_clean.ID,all_Subds(nS));
    
    data_clean.PowDelta(this_line)=all_Pow(nS,2);
    data_clean.PowTheta(this_line)=all_Pow(nS,3);
    data_clean.PowAlpha(this_line)=all_Pow(nS,4);
    
    data_clean.PowDeltaN(this_line)=all_Pow(nS,2)-all_Pow(nS,5);
    data_clean.PowThetaN(this_line)=all_Pow(nS,3)-all_Pow(nS,5);
    data_clean.PowAlphaN(this_line)=all_Pow(nS,4)-all_Pow(nS,5);
    
    data_clean.PowDeltaN_W(this_line)=all_Pow(nS,12)-all_Pow(nS,15);
    data_clean.PowThetaN_W(this_line)=all_Pow(nS,13)-all_Pow(nS,15);
    data_clean.PowAlphaN_W(this_line)=all_Pow(nS,14)-all_Pow(nS,15);
    
    data_clean.PowThetaAlpha(this_line)=all_Pow(nS,3)-all_Pow(nS,4);
    
    data_clean.PowSlope(this_line)=all_Pow(nS,7);
    data_clean.PowBG(this_line)=all_Pow(nS,6);
    
end

writetable(data_clean,'../Edison_Tables/Clean_Data_SS_Pow.txt');

%% Slept or Not
figure;
temp=data_clean.InsightPost(data_clean.SleepEdison==0);
simpleBarPlot(1,temp,[1 1 1;0 0 0],0.85,'r');

temp=data_clean.InsightPost(data_clean.SleepEdison==1);
simpleBarPlot(2,temp,[.7 .7 .7;0 0 0],0.85,'r');
xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'No','Yes'});
xlabel('Slept')
ylabel('Prob Insight')
format_fig

%% Wake N1 or N2
figure;
temp=data_clean.InsightPost(data_clean.SleepEdison==0);
simpleBarPlot(1,temp,[1 1 1;0 0 0],0.85,'r');

temp=data_clean.InsightPost(data_clean.SleepEdison==1 & data_clean.HypnoEdisonN2==0);
simpleBarPlot(2,temp,[.7 .7 .7;0 0 0],0.85,'r');

temp=data_clean.InsightPost(data_clean.SleepEdison==1 & data_clean.HypnoEdisonN2~=0);
simpleBarPlot(3,temp,[0 0 0;0 0 0],0.85,'r');


xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'No','N1','N2'});
xlabel('Slept')
ylabel('Prob Insight')

format_fig


%% Slept or Not - Delta Power
varname='PowDeltaN';
figure;
temp=data_clean.(varname)(data_clean.SleepEdison==0 & data_clean.InsightPost==0);
simpleBarPlot(1-0.2,temp,[1 1 1;1 0 0],0.35,'k');

temp=data_clean.(varname)(data_clean.SleepEdison==0 & data_clean.InsightPost==1);
simpleBarPlot(1+0.2,temp,[.7 .7 .7;1 0 0],0.35,'k');


temp=data_clean.(varname)(data_clean.SleepEdison==1 & data_clean.InsightPost==0);
simpleBarPlot(2-0.2,temp,[1 1 1;0 0 1],0.35,'k');

temp=data_clean.(varname)(data_clean.SleepEdison==1 & data_clean.InsightPost==1);
simpleBarPlot(2+0.2,temp,[.7 .7 .7;0 0 1],0.35,'k');

xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'No','Yes'});
xlabel('Slept')
ylabel('Delta Power')
format_fig

%% Slept or Not - Delta Power
figure;
temp=data_clean.PowDeltaN_W(data_clean.SleepEdison==0);
simpleBarPlot(1,temp,[1 1 1;0 0 0],0.85,'r');

temp=data_clean.PowDeltaN_W(data_clean.SleepEdison==1 & data_clean.HypnoEdisonN2==0);
simpleBarPlot(2,temp,[.7 .7 .7;0 0 0],0.85,'r');

temp=data_clean.PowDeltaN_W(data_clean.SleepEdison==1 & data_clean.HypnoEdisonN2~=0);
simpleBarPlot(3,temp,[0 0 0;0 0 0],0.85,'r');


xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'No','N1','N2'});
xlabel('Slept')
ylabel('Delta Power')

format_fig