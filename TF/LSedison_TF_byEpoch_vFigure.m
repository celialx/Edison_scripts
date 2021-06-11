%%
clear all
close all

% Fieldtrip
path_fieldtrip='/Users/thandrillon/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

% Thomas' toolkit

% FOOOF
fooof_path='/Users/thandrillon/WorkGit/projects/ext/fooof_mat/';
f_range = [2, 30];
settings = struct();  % Use defaults
addpath(genpath(fooof_path));

path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

data_path='/Users/thandrillon/Data/LS_Edison/EDF_fixed';
files=dir([data_path filesep '*.edf']);

<<<<<<< HEAD
path_export='/Users/thandrillon/Work/local/export_fig/';
=======
path_export='/Users/tand0009/Work/local/export_fig/';
>>>>>>> 5757eba8769d25b9c6d48a73c88db32511e24656
addpath(genpath(path_export));
ColorsGroup=[93 175 117;
    133 189 181;
    67 115 128]/256;

ColorsSolver=[196 146 46 ; 177 68 22]/256;
% % Wake = 55 179 111
% % N1 = 115 191 181
% % N2 = 50 116 130
% % Solvers = 206 144 45
% % Nonsolvers = 192 50 2


%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).

%%
nc=0;
StageDrop=[];
all_Pow=[];
all_Subds=[];
all_ERP=[];
all_AlphaPeaks=[];

all_Power_Drops=[];
all_Subds_Drops=[];

ParamExport_SubID=cell(0,0);
ParamExport=[];

all_maxAbs=[];

for nF=1:length(files)
    %%% Select subject and file
    File_Name       = files(nF).name;
    Folder_Name     = files(nF).folder;
    SubdID          = File_Name;
    sep             = findstr(File_Name,'_');
    SubdID          = SubdID(1:sep(1)-1);
    
%     if ~strcmp(SubdID,'03JF')
%         continue;
%     end
    if exist([Folder_Name filesep 'T_' SubdID '.mat'])==0
        warning('matrix file missing');
        continue;
    end
    load([Folder_Name filesep 'T_' SubdID '.mat'])
    fprintf('... processing %s\n',File_Name);
    
    %%% Import data, select for breaks and cut for 30s epochs
    cfg=[];
    cfg.trialfun        = 'LSedison_trialfun_break_contepochs';
    cfg.dataset         = [Folder_Name filesep  File_Name];
    cfg.demean          = 'yes';
    cfg.T               = T;
    cfg.lengthepochs    = 30;
    cfg                 = ft_definetrial(cfg);
    cfg.channel         = {'Fp1-A2','C3-A2','O1-A2'};
    data                = ft_preprocessing(cfg); % read raw data
    hdr                 = ft_read_header([Folder_Name filesep File_Name]);
    
    ParamExport=[ParamExport ; hdr.orig.PhysMin(1:3) hdr.orig.PhysMax(1:3)];
    ParamExport_SubID=[ParamExport_SubID ; {SubdID}];
    %%% Reject trials with trials above a max amplitude threshold of 150uV
    maxAbs = [];
    for k=1:length(data.trial)
        maxAbs(k,:)     = max(abs(data.trial{k}),[],2);
    end
    all_maxAbs=[all_maxAbs ; maxAbs];
    threshArt           = 150;
    if (100*mean(maxAbs(:,3)>threshArt))>50
        warning('Too much trials above Thr. Will Skip');
        continue;
    end
    nc=nc+1;
    pption_rejTr(nc,:)  = [sum(maxAbs(:,3)>threshArt) size(maxAbs,1) 100*mean(maxAbs(:,3)>threshArt)];
    fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(maxAbs(:,3)>threshArt),100*mean(maxAbs(:,3)>threshArt),threshArt)
    
    %%% Extract power using the Welch's technique
    TFRhann=[];
    for nCh=1:length(data.label)
        for nTr=1:length(data.trial)
            w_window                        = 6*data.fsample;
            w_overlap                       = w_window/2;
            df                              = 0.2;
            freqV                           = 1:0.2:30;
            signal                          = data.trial{nTr}(nCh,:);
            [pow,faxis]                     = pwelch(signal,w_window,w_overlap,freqV,data.fsample,'psd');
            TFRhann.powspctrm(nTr,nCh,:,1)  = pow;
            TFRhann.signal(nTr,nCh,:,1)     = signal;
        end
    end
    freqs = faxis;

    %%%% Extract power 5 minutes before drop
    AASM_Score=T.Stade;
    Beg_Task=(T.Epoch(find(T.Start(:,1))));
    End_Task=T.Epoch(find(T.End(:,1)));
    Drops=find(T.Ball(:,1));
    DropEpochs=T.Epoch(find(T.Ball(:,1)));
    if ~isempty(Drops)
        cat_data=[];
        DropsRel=Drops(1)-find(T.Start(:,1));
        for nTr=1:length(data.trial)
            cat_data=[cat_data  data.trial{nTr}(3,:)];
        end
        if DropsRel*hdr.Fs-5*60*hdr.Fs>1 & DropsRel*hdr.Fs+60*hdr.Fs<length(cat_data)
            cut_data=cat_data(DropsRel*hdr.Fs-5*60*hdr.Fs:DropsRel*hdr.Fs+60*hdr.Fs);
            
            w_window=6*data.fsample;
            w_overlap=w_window/2;
            df=0.2;
            freqV=1:0.2:30;
            nwindow=length(cut_data)/w_window;
            m=1;
            cut_pow=[];
            while (m+12*hdr.Fs)<length(cut_data)
                signal=cut_data(m:m+12*hdr.Fs);
                [pow,faxis_cut] = pwelch(signal,w_window,w_overlap,freqV,data.fsample,'psd');
                m=m+3*hdr.Fs;
                cut_pow(end+1,:)=10*log10(pow)';
            end
            all_Power_Drops=cat(3,all_Power_Drops,cut_pow);
            all_Subds_Drops=[all_Subds_Drops ; {SubdID}];
        end
    end
    
    %%% Match 30s epochs with scoring
%     % Uncomment if you want to focus on the epochs before the first drop
%     if ~isempty(Drops)
%         AASM_Score(Drops(1):end)=NaN;
%     end
    AASM_Score          = AASM_Score(T.Epoch>=Beg_Task & T.Epoch<=End_Task);
    AASM_Score_Epochs   = (T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score_Epochs2  = unique(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score2         = nan(1,length(data.trial));
    for j=1:length(data.trial)
        if sum(isnan(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)))))==length(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j))))
            AASM_Score2(j) = nan;
        else
            temp=unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)));
            AASM_Score2(j) = temp(1);
        end
    end
    % mark epochs to be rejected
    AASM_Score2(find(maxAbs(:,3)>threshArt)) =  NaN;
    
    
    %%% Apply FOOOF on spectrum and gather data
    Stages=0:2; % keep only wake, N1 and N2 stages (no NaN stage)
    fooof_results_perStage=[];
    for nSta=1:4
        if nSta<4
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
            if sum(isnan(temp))==length(squeeze(temp)')
                fooof_results_perStage{nSta}.background_params = nan(1,2);
                
                FOOOF_Pow_perStage(nc,nSta,:)   = nan(1,141);
                FOOOF2_Pow_perStage(nc,nSta,:)  = nan(1,141);
                
            else
                fooof_results_perStage{nSta}    = fooof(freqs,temp', f_range, settings,1);
                
                FOOOF_Pow_perStage(nc,nSta,:)   = fooof_results_perStage{nSta}.fooofed_spectrum;
                FOOOF2_Pow_perStage(nc,nSta,:)  = fooof_results_perStage{nSta}.fooofed_spectrum- fooof_results_perStage{nSta}.ap_fit;
            end
            Pow_perStage(nc,nSta,:)             = squeeze(nanmean(mean(log(TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
        else % gather all stages
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(ismember(AASM_Score2,Stages(2:end)),2,:,:)),1),4));
            if sum(isnan(temp))==length(squeeze(temp)')
                fooof_results_perStage{nSta}.background_params = nan(1,2);
                
                FOOOF_Pow_perStage(nc,nSta,:)   = nan(1,141);
                FOOOF2_Pow_perStage(nc,nSta,:)  = nan(1,141);
            else
                fooof_results_perStage{nSta}    = fooof(freqs,temp', f_range, settings,1);
                FOOOF_Pow_perStage(nc,nSta,:)   = fooof_results_perStage{nSta}.fooofed_spectrum;
                FOOOF2_Pow_perStage(nc,nSta,:)  = fooof_results_perStage{nSta}.fooofed_spectrum- fooof_results_perStage{nSta}.ap_fit;
            end
            Pow_perStage(nc,nSta,:)             = squeeze(nanmean(mean(log(TFRhann.powspctrm(ismember(AASM_Score2,Stages(2:end)),2,:,:)),1),4));
        end
    end
    
    if ~isempty(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')))
        Pow_all(nc,:)                 = squeeze(nanmean(nanmean(log(TFRhann.powspctrm(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')),3,:,:)),1),4));
        fooof_results                 = fooof(freqs,squeeze(nanmean(nanmean((TFRhann.powspctrm(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')),3,:,:)),1),4))', f_range, settings,1);
        
        alpha_peaks         = fooof_results.peak_params(fooof_results.peak_params(:,1)>7 & fooof_results.peak_params(:,1)<12,:);
        if ~isempty(alpha_peaks)
            alpha_peaks     = alpha_peaks(alpha_peaks(:,2)==max(alpha_peaks(:,2)),:);
            all_AlphaPeaks  = [all_AlphaPeaks ; [nF sum(AASM_Score==0) sum(AASM_Score==1) sum(AASM_Score==2) alpha_peaks]];
        else
            all_AlphaPeaks  = [all_AlphaPeaks ; [nF sum(AASM_Score==0) sum(AASM_Score==1) sum(AASM_Score==2) nan(1,3)]];
        end
    else
        Pow_all(nc,:)                 = nan;
        fooof_results.peak_params       = [];
        fooof_results.background_params = nan(1,2);
    end
    for nCh=1:3
        if ~isempty(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')))
            Pow_allE(nc,nCh,:)        = squeeze(nanmean(nanmean(log(TFRhann.powspctrm(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')),nCh,:,:)),1),4));
            fooof_results2            = fooof(freqs,squeeze(nanmean(nanmean((TFRhann.powspctrm(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')),nCh,:,:)),1),4))', f_range, settings,1);
            FOOOF_Pow_allE(nc,nCh,:)  = fooof_results2.fooofed_spectrum;
            FOOOF2_Pow_allE(nc,nCh,:) = fooof_results2.fooofed_spectrum-fooof_results2.ap_fit;
        else
            Pow_allE(nc,nCh,:)        = nan;
            FOOOF_Pow_allE(nc,nCh,:)  = nan;;
            FOOOF2_Pow_allE(nc,nCh,:) = nan;;
        end
    end
    
    all_Pow = [all_Pow ; [nF mean(Pow_all(nc,freqs>2 & freqs<4)) mean(Pow_all(nc,freqs>4 & freqs<8)) mean(Pow_all(nc,freqs>8 & freqs<11)) mean(Pow_all(nc,freqs>1 & freqs<40)) fooof_results.background_params length(AASM_Score2) sum(AASM_Score2==0) sum(AASM_Score2==1) sum(AASM_Score2==2) ...
        mean((Pow_perStage(nc,1,freqs>2 & freqs<4))) mean((Pow_perStage(nc,1,freqs>4 & freqs<8))) mean((Pow_perStage(nc,1,freqs>8 & freqs<11))) mean((Pow_perStage(nc,1,freqs>1 & freqs<40))) fooof_results_perStage{1}.background_params...
        mean((Pow_perStage(nc,2,freqs>2 & freqs<4))) mean((Pow_perStage(nc,2,freqs>4 & freqs<8))) mean((Pow_perStage(nc,2,freqs>8 & freqs<11))) mean((Pow_perStage(nc,2,freqs>1 & freqs<40))) fooof_results_perStage{2}.background_params...
        mean((Pow_perStage(nc,3,freqs>2 & freqs<4))) mean((Pow_perStage(nc,3,freqs>4 & freqs<8))) mean((Pow_perStage(nc,3,freqs>8 & freqs<11))) mean((Pow_perStage(nc,3,freqs>1 & freqs<40))) fooof_results_perStage{3}.background_params...
        mean((Pow_perStage(nc,4,freqs>2 & freqs<4))) mean((Pow_perStage(nc,4,freqs>4 & freqs<8))) mean((Pow_perStage(nc,4,freqs>8 & freqs<11))) mean((Pow_perStage(nc,4,freqs>1 & freqs<40))) fooof_results_perStage{4}.background_params ]];
    all_Subds = [all_Subds ; {SubdID}];
end
freqs2=fooof_results.freqs;

%%
Freqs=[3.2 4.4; 7 9; 9 9.8];
data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');
data_clean.PowDelta=nan(size(data_clean,1),1);
data_clean.PowTheta=nan(size(data_clean,1),1);
data_clean.PowAlpha=nan(size(data_clean,1),1);

data_clean.PowDelta_dev=nan(size(data_clean,1),1);
data_clean.PowTheta_dev=nan(size(data_clean,1),1);
data_clean.PowAlpha_dev=nan(size(data_clean,1),1);

data_clean.PowThetaAlpha=nan(size(data_clean,1),1);
data_clean.PowThetaDelta=nan(size(data_clean,1),1);


data_clean.PowDeltaN=nan(size(data_clean,1),1);
data_clean.PowThetaN=nan(size(data_clean,1),1);
data_clean.PowAlphaN=nan(size(data_clean,1),1);

data_clean.PowDeltaN_W=nan(size(data_clean,1),1);
data_clean.PowThetaN_W=nan(size(data_clean,1),1);
data_clean.PowAlphaN_W=nan(size(data_clean,1),1);

data_clean.PowDeltaN_S=nan(size(data_clean,1),1);
data_clean.PowThetaN_S=nan(size(data_clean,1),1);
data_clean.PowAlphaN_S=nan(size(data_clean,1),1);

data_clean.AlphaPeakFreq=nan(size(data_clean,1),1);
data_clean.AlphaPeakAmp=nan(size(data_clean,1),1);

data_clean.PowSlope=nan(size(data_clean,1),1);
data_clean.PowBG=nan(size(data_clean,1),1);

data_clean.SleepGroup=nan(size(data_clean,1),1);
data_clean.SleepGroup(data_clean.HypnoEdisonN1==0  & data_clean.HypnoEdisonN2==0)=0;
data_clean.SleepGroup(data_clean.HypnoEdisonN1~=0  & data_clean.HypnoEdisonN2==0)=1;
data_clean.SleepGroup(data_clean.HypnoEdisonN1~=0  & data_clean.HypnoEdisonN2~=0)=2;
data_clean.SleepGroup=categorical(data_clean.SleepGroup);
% data_clean.SleepGroup=reordercats(data_clean.SleepGroup,[2 1 3]);


All_Conds=[];
All_Insight=[];
for nS=1:length(all_Subds)
    this_line=match_str(data_clean.ID,all_Subds(nS));
    
    data_clean.PowDelta(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>=Freqs(1,1) & freqs2<=Freqs(1,2)),3));
    data_clean.PowTheta(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>=Freqs(2,1) & freqs2<=Freqs(2,2)),3));
    data_clean.PowAlpha(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>=Freqs(3,1) & freqs2<=Freqs(3,2)),3));
    
    data_clean.PowDelta_dev(this_line)=squeeze(mean(FOOOF2_Pow_allE(nS,3,freqs2>=Freqs(1,1) & freqs2<=Freqs(1,2)),3));
    data_clean.PowTheta_dev(this_line)=squeeze(mean(FOOOF2_Pow_allE(nS,3,freqs2>=Freqs(2,1) & freqs2<=Freqs(2,2)),3));
    data_clean.PowAlpha_dev(this_line)=squeeze(mean(FOOOF2_Pow_allE(nS,3,freqs2>=Freqs(3,1) & freqs2<=Freqs(3,2)),3));
    
    
    data_clean.PowDeltaN(this_line)=all_Pow(nS,2)-all_Pow(nS,5);
    data_clean.PowThetaN(this_line)=all_Pow(nS,3)-all_Pow(nS,5);
    data_clean.PowAlphaN(this_line)=all_Pow(nS,4)-all_Pow(nS,5);
    
    data_clean.PowDeltaN_W(this_line)=all_Pow(nS,12)-all_Pow(nS,15);
    data_clean.PowThetaN_W(this_line)=all_Pow(nS,13)-all_Pow(nS,15);
    data_clean.PowAlphaN_W(this_line)=all_Pow(nS,14)-all_Pow(nS,15);
    
    data_clean.PowDeltaN_S(this_line)=all_Pow(nS,30)-all_Pow(nS,33);
    data_clean.PowThetaN_S(this_line)=all_Pow(nS,31)-all_Pow(nS,33);
    data_clean.PowAlphaN_S(this_line)=all_Pow(nS,32)-all_Pow(nS,33);
    
    data_clean.PowSlope(this_line)=all_Pow(nS,7);
    data_clean.PowBG(this_line)=all_Pow(nS,6);
    
    data_clean.AlphaPeakFreq(this_line)=all_AlphaPeaks(nS,5);
    data_clean.AlphaPeakAmp(this_line)=all_AlphaPeaks(nS,6);
    
    if isempty(double(data_clean.SleepGroup(this_line)))
        All_Conds=[All_Conds ;  NaN];
        All_Insight=[All_Insight ;  NaN];
    else
        All_Conds=[All_Conds ;  double(data_clean.SleepGroup(this_line))];
        All_Insight=[All_Insight ;  (data_clean.InsightPost(this_line))];
    end
end
All_Conds2=double(All_Conds~=1);
All_Conds2(find(isnan(All_Conds)))=nan;
% writetable(data_clean,'../Edison_Tables/Clean_Data_SS_Pow.txt');


%%
figure; hold on;
set(gcf,'Position',[440   255   338   543]);
for nCond=1:3
    temp_toplot=squeeze(Pow_perStage(:,nCond,:));
    %     temp_toplot=(temp_toplot)-repmat(nanmean(temp_toplot,2),1,length(freqs));
    plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
end
xlabel('Frequency (Hz)')
legend({'WK','N1','N2'})
format_fig;
xlim([2 20])
title({'Power Spectrum','per Sleep Stage'})
% ylim([-1 3])

%%%
figure; hold on;
set(gcf,'Position',[779   253   338   545]);
myConds={[1],[2 3]};
freqs=freqs;
for nCond=1:3
    temp_toplot=squeeze(Pow_allE(All_Conds==nCond,3,:));
    %     temp_toplot=(temp_toplot)-repmat(nanmean(temp_toplot,2),1,length(freqs));
    plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
end
xlabel('Frequency (Hz)')
legend({'WK','N1','N2'})
format_fig;
xlim([2 20])
title({'Power Spectrum','per Group'})
% ylim([-1 3])

%%
% [realpos realneg]=get_cluster_permutation_aov(squeeze(Pow_allE((~isnan(All_Conds)),3,:)),[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],...
%     0.05,0.1,100,freqs);
redo=0;
if redo
    temp_power=squeeze(FOOOF_Pow_allE((~isnan(All_Conds)),3,:));
    [realpos_lin realneg_lin]=get_cluster_permutation_aov(zscore(temp_power),[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],...
        0.05,0.1,1000,freqs2);
    [realpos_quad realneg_quad]=get_cluster_permutation_aov(sqrt(zscore(temp_power).^2),[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],...
        0.05,0.1,1000,freqs2);
    save('result_clusterperm_onfreq_ANOVA_FOOOFpower','realpos_lin','realneg_lin','realpos_quad','realneg_quad')
else
    load('result_clusterperm_onfreq_ANOVA_FOOOFpower')
end
%%
    temp_power=squeeze(FOOOF_Pow_allE((~isnan(All_Conds)),3,:));
    myfreqs=freqs2;
figure;
subplot(1,3,1);
hold on;
for nCond=1:2
    for nCond2=1:2
        temp_toplot=squeeze(temp_power(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Insight(~isnan(All_Conds))==(nCond2-1),:));
        if nCond2==1
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 20])
% ylim([-1 4])

legend({'WK | NO', 'WK | YS','SL | NO','SL | YE','sleep','insight','interaction'})

subplot(1,3,2);
temp_power_zscore=zscore(temp_power);
hold on;
for nCond=1:2
    for nCond2=1:2
        temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Insight(~isnan(All_Conds))==(nCond2-1),:);
        if nCond2==1
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 20])
% ylim([-1 4])

hold on;
scatter(myfreqs(find(realpos_lin{1}.clusters)),-1+0.2*ones(1,length(find(realpos_lin{1}.clusters))),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(myfreqs(find(realpos_lin{2}.clusters)),-1.5+0.2*ones(1,length(find(realpos_lin{2}.clusters))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(myfreqs(find(realpos_lin{3}.clusters)),-2+0.2*ones(1,length(find(realpos_lin{3}.clusters))),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');


subplot(1,3,3);
temp_power_zscore=zscore(temp_power).^2;
hold on;
for nCond=1:2
    for nCond2=1:2
        temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Insight(~isnan(All_Conds))==(nCond2-1),:);
        if nCond2==1
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 20])
% ylim([-1 4])

hold on;
scatter(myfreqs(find(realpos_quad{1}.clusters)),-1+0.2*ones(1,length(find(realpos_quad{1}.clusters))),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(myfreqs(find(realpos_quad{2}.clusters)),-1.5+0.2*ones(1,length(find(realpos_quad{2}.clusters))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(myfreqs(find(realpos_quad{3}.clusters)),-2+0.2*ones(1,length(find(realpos_quad{3}.clusters))),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');


%% Box plot Linear effects
%%%% clean clusters
realpos_quad2=realpos_quad;
realneg_quad2=realneg_quad;
realpos_lin2=realpos_lin;
realneg_lin2=realneg_lin;
for k=1:length(realpos_quad2)
    realpos_quad2{k}.clusters2=realpos_quad2{k}.clusters;
    for m=unique(realpos_quad2{k}.clusters2)
        if sum(realpos_quad2{k}.clusters2==m)<5
            realpos_quad2{k}.clusters2(realpos_quad2{k}.clusters2==m)=0;
        end
    end
    realneg_quad2{k}.clusters2=realneg_quad2{k}.clusters;
    for m=unique(realneg_quad2{k}.clusters2)
        if sum(realneg_quad2{k}.clusters2==m)<5
            realneg_quad2{k}.clusters2(realneg_quad2{k}.clusters2==m)=0;
        end
    end
    realpos_lin2{k}.clusters2=realpos_lin2{k}.clusters;
    for m=unique(realpos_lin2{k}.clusters2)
        if sum(realpos_lin2{k}.clusters2==m)<5
            realpos_lin2{k}.clusters2(realpos_lin2{k}.clusters2==m)=0;
        end
    end
    realneg_lin2{k}.clusters2=realneg_lin2{k}.clusters;
    for m=unique(realneg_lin2{k}.clusters2)
        if sum(realneg_lin2{k}.clusters2==m)<5
            realneg_lin2{k}.clusters2(realneg_lin2{k}.clusters2==m)=0;
        end
    end
end

Freqs=[3.2 4.4; 9 9.8];
figure;
hold on;
cmap=cbrewer('seq','Blues',16); 
ylims=[-0.4 1.2];
for nFreq=1:size(Freqs,1)
    patch([Freqs(nFreq,1) Freqs(nFreq,2) Freqs(nFreq,2) Freqs(nFreq,1) Freqs(nFreq,1)],[ylims(1) ylims(1) ylims(2) ylims(2) ylims(1)],...
        cmap(nFreq*4,:),'EdgeColor',cmap(nFreq*4,:),'FaceAlpha',0.5);
end

ColorsGroup2=ColorsGroup([1 3],:);
set(gcf,'Position',[1     1   460   804]);
hp=[];
for nCond=1:2
    for nCond2=1:2
        temp_toplot=squeeze(FOOOF_Pow_allE(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),3,:));
        if nCond2==1
            hp(end+1)=plot(freqs2,nanmean(temp_toplot),'Color',ColorsGroup2(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            hp(end+1)=plot(freqs2,nanmean(temp_toplot),'Color',ColorsGroup2(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([2 20])
ylabel('Power (dB)')
% 
set(gca,'LineWidth',2);

scatter(myfreqs(find(realpos_lin2{2}.clusters2)),-0.22*ones(1,length(find(realpos_lin2{2}.clusters2))),'Marker','s','MarkerEdgeColor',[1 0 0]*0.9,'MarkerFaceColor',[1 0 0]*0.9,'SizeData',64);
scatter(myfreqs(find(realpos_quad2{2}.clusters2)),-0.2*ones(1,length(find(realpos_quad2{2}.clusters2))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','SizeData',64);

% scatter(myfreqs(find(realpos_lin{1}.clusters)), 1.08*ones(1,length(find(realpos_lin{1}.clusters))),'Marker','o','MarkerEdgeColor',[1 1 1]*0.2,'MarkerFaceColor',[1 1 1]*0.2);
% scatter(myfreqs(find(realpos_quad{1}.clusters)),1.10*ones(1,length(find(realpos_quad{1}.clusters))),'Marker','o','MarkerEdgeColor',[1 1 1]*0.5,'MarkerFaceColor',[1 1 1]*0.5);
hl=legend(hp,{'WK | NO', 'WK | YS','SL | NO','SL | YE'},'Position',[0.7174    0.6405    0.2293    0.1063]);

export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_TimePlot.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_TimePlot.png'],'-r 300')

% Freqs=[1 4; 4 8; 8 12 ; 15 20];
TitleFreqs={'Delta','Alpha'};
figure;
set(gcf,'Position',[462     1   348   804]);
for nFreq=1
    subplot(size(Freqs,1),1,1);
    format_fig;
    for nCond=1:length(myConds)
        for nCond2=1:2
            temp_toplot=squeeze(FOOOF_Pow_allE(:,3,:));
            %         temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
            temp_toplot=mean(temp_toplot(:,freqs2>=Freqs(nFreq,1) & freqs2<=Freqs(nFreq,2)),2);
            temp_toplot=(nanzscore(temp_toplot));
            temp_toplot=temp_toplot(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),:);
            %%%%% Plot
            hold on;
            Pos=2*nCond2+0.2*(2*nCond-3); data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
            if nCond2==1
                colorBar=[ColorsGroup2(nCond,:) ;ColorsGroup2(nCond,:)];
            else
                colorBar=[ColorsGroup2(nCond,:) ; 0 0 0];
            end
            line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
            line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
            line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
            line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
            line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
            
            patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
                [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
            
            xspread=(rand(1,length(data))-0.5)*widthBar/4+2*nCond2+.65*(2*nCond-3);
            scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
        end
    end
    xlim([.5 5.5])
    set(gca,'XTick',[2 4],'XTickLabel',{'NO','YES'},'LineWidth',2)
    ylabel(sprintf('z(Power %s)',TitleFreqs{nFreq}))
    xlabel('Insight')
end

%%% Box plot Quadratic effects
for nFreq=2
     subplot(size(Freqs,1),1,2);
   format_fig;
    for nCond=1:length(myConds)
        for nCond2=1:2
            temp_toplot=squeeze(Pow_allE(:,3,:));
            %         temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
            temp_toplot=mean(temp_toplot(:,freqs>=Freqs(nFreq,1) & freqs<=Freqs(nFreq,2)),2);
            temp_toplot=sqrt(nanzscore(temp_toplot).^2);
            temp_toplot=temp_toplot(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),:);
            %%%%% Plot
            hold on;
            Pos=2*nCond2+0.2*(2*nCond-3); data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
            if nCond2==1
                colorBar=[ColorsGroup2(nCond,:) ;ColorsGroup2(nCond,:)];
            else
                colorBar=[ColorsGroup2(nCond,:) ; 0 0 0];
            end
            line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
            line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
            line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
            line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
            line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
            
            patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
                [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
            
            xspread=(rand(1,length(data))-0.5)*widthBar/4+2*nCond2+.65*(2*nCond-3);
            scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
        end
    end
    xlim([.5 5.5])
    set(gca,'XTick',[2 4],'XTickLabel',{'NO','YES'},'LineWidth',2)
    ylabel(sprintf('z(Power %s)^2',TitleFreqs{nFreq}))
    xlabel('Insight')
end
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_DotPlot.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_DotPlot.png'],'-r 300')

%%
figure;
stepbin=33;
clear bin_*
titlePlots={'Delta','Alpha'};
set(gcf,'Position',[462     1   348   804]);
for nplot=1:2
    subplot(2,1,nplot);
    hold on;
    format_fig;
    if nplot==1
        tempX=nanzscore(data_clean.PowDelta);
    elseif nplot==2
        tempX=nanzscore(data_clean.PowAlpha);
    end
%         tempY=data_clean.RT_fromAha-data_clean.RT_beforeAha;
    tempY=data_clean.InsightPost;
    tempC=double(data_clean.SleepGroup~='0')+1;
    tempX(isnan(tempC))=[];
    tempY(isnan(tempC))=[];
    tempC(isnan(tempC))=[];
    
    bins=prctile(tempX,0:stepbin:100);
    bin_values=[];
    bin_values{1,1}=(tempX(tempX<bins(2)));
    bin_values{1,2}=(tempY(tempX<bins(2)));
    bin_values{1,3}=(tempC(tempX<bins(2)));
    for k=2:length(bins)-2
        bin_values{k,1}=(tempX(tempX>=bins(k) & tempX<=bins(k+1)));
        bin_values{k,2}=(tempY(tempX>=bins(k) & tempX<=bins(k+1)));
        bin_values{k,3}=(tempC(tempX>=bins(k) & tempX<=bins(k+1)));
    end
    bin_values{length(bins)-1,1}=(tempX(tempX>bins(length(bins)-1)));
    bin_values{length(bins)-1,2}=(tempY(tempX>bins(length(bins)-1)));
    bin_values{length(bins)-1,3}=(tempC(tempX>bins(length(bins)-1)));
    
    bin_values2=[];
    bin_values_sem=[];
    for kbin=1:size(bin_values,1)
        Pos=kbin; data=bin_values{kbin,2}; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
        %         colorBar=[.5 .5 .5; 0 0 0];
        colorBar=[];
        for l=1:length(bin_values{kbin,3})
            if isnan(bin_values{kbin,3}(l))
                colorBar(l,:)=[0.5 0.5 0.5];
            else
                colorBar(l,:)=ColorsGroup2(bin_values{kbin,3}(l),:);
            end
        end
        xspread=(rand(1,length(data))-0.5)*widthBar/5+Pos;
        yspread=(rand(1,length(data))-0.5)*widthBar/10+data';
        scatter(xspread',yspread',[],colorBar,'filled','Marker',markerType,'MarkerFaceAlpha',0.75,'SizeData',sizeDot/3);
        
        bin_values2(kbin,1)=mean(data);
        bin_values_sem(kbin,1)=sem(data);
    end
    errorbar(1:length(bins)-1,bin_values2(:,1),-bin_values_sem(:,1),+bin_values_sem(:,1),'Color','k','LineWidth',3);
    scatter(1:length(bins)-1,bin_values2(:,1),'Marker','o','SizeData',244,'MarkerFaceColor',[1 1 1]*0.7,'MarkerEdgeAlpha',.7,'MarkerEdgeColor','k','LineWidth',3);
    
    %     ylim([0 0.8])
    xlim([0.5 length(bins)-1+0.5])
    ylim([-0.08 1.08])
<<<<<<< HEAD
    set(gca,'XTick',1:length(bins),'XColor','k','YColor','k'); %,'XTickLabel',{'low','med','high'});
=======
    set(gca,'XTick',1:length(bins),'XColor','k','YColor','k','LineWidth',2); %,'XTickLabel',{'low','med','high'});
>>>>>>> 5757eba8769d25b9c6d48a73c88db32511e24656
    %     title(titlePlots{nplot})
    ylabel('Insight','Color','k')
    xlabel(['Power Bin ',titlePlots{nplot}],'Color','k')
end
%


export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_BinPlot.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_BinPlot.png'],'-r 300')


%%
anovan(nanzscore(data_clean.PowDelta),[data_clean.SleepEdison data_clean.InsightPost],'varnames',{'Sleep','Insight'},'model','full');
anovan(nanzscore(data_clean.PowAlpha).^2,[data_clean.SleepEdison data_clean.InsightPost],'varnames',{'Sleep','Insight'},'model','full');

%% Before Drops
%%% Match files
% all_Power_Drops=cat(3,all_Power_Drops,cut_pow);
% all_Subds_Drops=[all_Subds_Drops ; {SubdID}];

All_Conds_Drops=[];
All_Insigh_Drops=[];
for nS=1:length(all_Subds_Drops)
    this_line=match_str(data_clean.ID,all_Subds_Drops(nS));
    
    if isempty(double(data_clean.SleepEdison(this_line)))
        All_Conds_Drops=[All_Conds_Drops ;  NaN];
        All_Insigh_Drops=[All_Insigh_Drops ;  NaN];
    else
        All_Conds_Drops=[All_Conds_Drops ;  double(data_clean.SleepEdison(this_line))];
        All_Insigh_Drops=[All_Insigh_Drops ;  (data_clean.InsightPost(this_line))];
    end
end

figure; 
set(gcf,'position',[440     2   577   796]);

Power_byTime=cell(1,size(Freqs,1));
xTime=-5*60:3:1*60-13;
nFc=0;
hp=[];
for nFreqs=[1:2]
    nFc=nFc+1;
    subplot(2,1,nFc);
    format_fig;
    
    for nTime=1:size(all_Power_Drops,1)
        temp_power=squeeze(nanmean(all_Power_Drops(nTime,faxis_cut>Freqs(nFreqs,1) & faxis_cut<Freqs(nFreqs,2),:),2));
        Power_byTime{nFreqs}(nTime,:)=temp_power;
    end
    
    temp_plot=Power_byTime{nFreqs}(:,All_Insigh_Drops==0);
[~, hp(1)]=simpleTplot(xTime,temp_plot',0,ColorsSolver(2,:),0,'-',0.5,1,4,1,2);
hold on;
temp_plot=Power_byTime{nFreqs}(:,All_Insigh_Drops==1);
[~, hp(2)]=simpleTplot(xTime,temp_plot',0,ColorsSolver(1,:),0,'-',0.5,1,4,1,2);

if nFreqs==1
    hl=legend(hp,{'Non-Solvers','Solvers'},'position',[0.1724    0.8606    0.2366    0.0559]);
end
ylabel({'Power',titlePlots{nFreqs}})
xlabel('Time from Drop (s)')
end
% 
% 
% % Est_byTime=cell(1,length(nFreqs));
% % for nTime=1:size(all_Power_Drops,1)
% %     table_cut=[All_Conds_Drops All_Insigh_Drops];
% %     temp_power=squeeze(mean(all_Power_Drops(nTime,faxis_cut>8 & faxis_cut<12,:),2)); %./squeeze(mean(all_Power_Drops(nTime,faxis_cut>3 & faxis_cut<7,:),2));
% %     table_cut=[table_cut temp_power];
% %     table_cut=array2table(table_cut,'VariableNames',{'SleepEdison','InsightPost','PowDeltaTheta'});
% %     mdl = fitglm(table_cut,'InsightPost~1+PowDeltaTheta','Distribution','binomial');
% %     for nFreqs=1:length(Freqs)
% %         Est_byTime(nTime,:,1)=table2array(mdl.Coefficients(:,1));
% %         Est_byTime(nTime,:,2)=table2array(mdl.Coefficients(:,2));
% %     end
% % end
% 
% figure;
% subplot(1,3,1); hold on;
% errorbar(-5*60:3:1*60-13,squeeze(Est_byTime(:,2,1)),squeeze(Est_byTime(:,2,2)));
% % hold on;
% % plot(squeeze(Est_byTime(:,3,1)));
