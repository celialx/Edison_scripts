%%
clear all;close all;

run localdef_spatialInsight.m

addpath(path_fieldtrip);
ft_defaults;

% % % FOOOF
% % fooof_path='C:\Users\OUDIETTE.Delphine\Documents\MATLAB\Toolboxs\fooof_mat-main\';
% % f_range = [2, 30];
% % settings = struct();  % Use defaults
% % addpath(genpath(fooof_path));

addpath(genpath(path_LSCPtools));

% natsort_path ='C:\Users\OUDIETTE.Delphine\Documents\MATLAB\Toolboxs\natsortfiles\';
% addpath(genpath(natsort_path));

files=dir([data_path filesep '*.edf']);

% To sort files by filename
% files = natsortfiles(files);

ColorsGroup=[55 179 111;
    115 191 181;
    50 116 130]/255;

ColorsSolver=[206 144 45 ; 192 50 2]/255;

%%
nc=0;
all_Pow=[];
all_Subds=[];

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
    
    
    if exist([Folder_Name filesep 'T_' SubdID '.mat'])==0
        warning('matrix file missing');
        continue;
    end
    load([Folder_Name filesep 'T_' SubdID '.mat'])
    fprintf('... processing %s\n',File_Name);
    
    %%% Import data, select for breaks and cut for 30s epochs
    cfg=[];
    cfg.trialfun        = 'LSinsight_trialfun_break_contepochs';
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
    % Fait uniquement sur la troisième électrode
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
    Beg_Task=(T.Epoch(find(T.Start(:,1)))); Beg_Task = Beg_Task(1);
    End_Task=T.Epoch(find(T.End(:,1))); End_Task = End_Task(1);
    
    
    %%% Match 30s epochs with scoring
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
    
    Stages=0:2; % keep only wake, N1 and N2 stages (no NaN stage)
    for nSta=1:4
        if nSta<4
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
            Pow_perStage(nc,nSta,:)             = squeeze(nanmean(mean(log(TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
        else % gather all stages
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(ismember(AASM_Score2,Stages(2:end)),2,:,:)),1),4));
            Pow_perStage(nc,nSta,:)             = squeeze(nanmean(mean(log(TFRhann.powspctrm(ismember(AASM_Score2,Stages(2:end)),2,:,:)),1),4));
        end
    end
    
    if ~isempty(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')))
        Pow_all(nc,:)                 = squeeze(nanmean(nanmean(log(TFRhann.powspctrm(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')),3,:,:)),1),4));
    else
        Pow_all(nc,:)                 = nan;
    end
    
    for nCh=1:3
        if ~isempty(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')))
            Pow_allE(nc,nCh,:)        = squeeze(nanmean(nanmean(log(TFRhann.powspctrm(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')),nCh,:,:)),1),4));
        else
            Pow_allE(nc,nCh,:)        = nan;
        end
    end
    
    %       all_Pow = [all_Pow ; [nF mean(Pow_all(nc,freqs>2 & freqs<4)) mean(Pow_all(nc,freqs>4 & freqs<8)) mean(Pow_all(nc,freqs>8 & freqs<11)) mean(Pow_all(nc,freqs>1 & freqs<40)) length(AASM_Score2) sum(AASM_Score2==0) sum(AASM_Score2==1) sum(AASM_Score2==2) ...
    %         mean((Pow_perStage(nc,1,freqs>2 & freqs<4))) mean((Pow_perStage(nc,1,freqs>4 & freqs<8))) mean((Pow_perStage(nc,1,freqs>8 & freqs<11))) mean((Pow_perStage(nc,1,freqs>1 & freqs<40)))...
    %         mean((Pow_perStage(nc,2,freqs>2 & freqs<4))) mean((Pow_perStage(nc,2,freqs>4 & freqs<8))) mean((Pow_perStage(nc,2,freqs>8 & freqs<11))) mean((Pow_perStage(nc,2,freqs>1 & freqs<40)))...
    %         mean((Pow_perStage(nc,3,freqs>2 & freqs<4))) mean((Pow_perStage(nc,3,freqs>4 & freqs<8))) mean((Pow_perStage(nc,3,freqs>8 & freqs<11))) mean((Pow_perStage(nc,3,freqs>1 & freqs<40)))...
    %         mean((Pow_perStage(nc,4,freqs>2 & freqs<4))) mean((Pow_perStage(nc,4,freqs>4 & freqs<8))) mean((Pow_perStage(nc,4,freqs>8 & freqs<11))) mean((Pow_perStage(nc,4,freqs>1 & freqs<40)))]];
    all_Subds = [all_Subds ; {SubdID}];
end

%%
Freqs=[1 3.6; 6.6 8.6; 20 30];

% Behav data;
load([Behav_Path 'All_Data_clean']);

Data.PowDelta=nan(size(Data,1),1);
Data.PowAlpha=nan(size(Data,1),1);
Data.PowBeta=nan(size(Data,1),1);

Data.SleepGroup=nan(size(Data,1),1);
Data.SleepGroup(Data.HypnoNapN1==0  & Data.HypnoNapN2==0)=0;
Data.SleepGroup(Data.HypnoNapN1~=0  & Data.HypnoNapN2==0)=1;
Data.SleepGroup(Data.HypnoNapN2~=0)=2;
Data.SleepGroup=categorical(Data.SleepGroup);

All_Conds=[];
All_Memory=[];
for nS=1:length(all_Subds)
    this_line=match_str(Data.Code,all_Subds(nS));
    if strcmp(Data.Code, all_Subds(nS))==0
        error('not the same ID');
    end
    Data.PowDelta(this_line)=squeeze(nanmean(Pow_allE(nS,3,freqs>=Freqs(1,1) & freqs<=Freqs(1,2)),3));
    Data.PowAlpha(this_line)=squeeze(nanmean(Pow_allE(nS,3,freqs>=Freqs(2,1) & freqs<=Freqs(2,2)),3));
    Data.PowBeta(this_line)=squeeze(nanmean(Pow_allE(nS,3,freqs>=Freqs(3,1) & freqs<=Freqs(3,2)),3));
    %     Data.RatioAT(this_line) = Data.PowAlpha(this_line)/Data.PowTheta(this_line);
    
    if isempty(double(Data.SleepGroup(this_line)))
        All_Conds=[All_Conds ;  NaN];
        All_Memory=[All_Memory ;  NaN];
    else
        All_Conds=[All_Conds ;  double(Data.SleepGroup(this_line))];
        All_Memory=[All_Memory ;  (Data.Corrprenotpost(this_line))];
    end
end

All_Conds2=double(All_Conds~=1);
All_Conds2(find(isnan(All_Conds)))=nan;

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
% [realpos realneg]=get_cluster_permutation_aov(squeeze(Pow_allE((~isnan(All_Conds)),3,:)),[All_Conds2(~isnan(All_Conds)) All_Memory(~isnan(All_Conds))],...
%     0.05,0.1,100,freqs);
redo=1;
if redo
    temp_power=squeeze(Pow_allE((~isnan(All_Conds)),3,:));
    group=[All_Conds2(~isnan(All_Conds2)) All_Memory(~isnan(All_Conds2))];
    totperm=1000;
    [realpos_lin realneg_lin]=get_cluster_permutation_lm(zscore(temp_power),group,{'Power','Cond','Mem'},'Power~1+Cond*Mem',0.05,0.1,totperm,faxis);
    [realpos_quad realneg_quad]=get_cluster_permutation_lm(zscore(temp_power).^2,group,{'Power','Cond','Mem'},'Power~1+Cond*Mem',0.05,0.1,totperm,faxis);

    save('result_clusterperm_Insight2_onfreq_LM','realpos_lin','realneg_lin','realpos_quad','realneg_quad', 'freqs')
else
    load('result_clusterperm_Insight2_onfreq_LM')
end
%%
temp_power=squeeze(Pow_allE((~isnan(All_Conds)),3,:));
myfreqs=freqs;
figure; set(gcf,'Position',[1         141        1000         656]);
subplot(1,3,1);
hold on;
for nCond=1:2
    for nCond2=1:2
        if nCond2==1
        temp_toplot=squeeze(temp_power(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Memory(~isnan(All_Conds))<median(All_Memory(~isnan(All_Conds))),:));
        else
        temp_toplot=squeeze(temp_power(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Memory(~isnan(All_Conds))>median(All_Memory(~isnan(All_Conds))),:));
        end
        if nCond2==1
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 30])
% ylim([-1 4])

legend({'WK | NO', 'WK | YS','SL | NO','SL | YE','sleep','insight','interaction'})

subplot(1,3,2);
temp_power_zscore=zscore(temp_power);
hold on;
for nCond=1:2
    for nCond2=1:2
        if nCond2==1
            temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Memory(~isnan(All_Conds))<median(All_Memory(~isnan(All_Conds))),:);
        else
            temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Memory(~isnan(All_Conds))>median(All_Memory(~isnan(All_Conds))),:);
        end
        if nCond2==1
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 30])
% ylim([-1 4])

hold on;
scatter(myfreqs(find(realpos_lin{1}.clusters)),-1+0.2*ones(1,length(find(realpos_lin{1}.clusters))),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(myfreqs(find(realpos_lin{2}.clusters)),-1.5+0.2*ones(1,length(find(realpos_lin{2}.clusters))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(myfreqs(find(realpos_lin{3}.clusters)),-2+0.2*ones(1,length(find(realpos_lin{3}.clusters))),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');

scatter(myfreqs(find(realneg_lin{1}.clusters)),-1.2+0.2*ones(1,length(find(realneg_lin{1}.clusters))),'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(myfreqs(find(realneg_lin{2}.clusters)),-1.7+0.2*ones(1,length(find(realneg_lin{2}.clusters))),'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(myfreqs(find(realneg_lin{3}.clusters)),-2+2.2*ones(1,length(find(realneg_lin{3}.clusters))),'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r');


subplot(1,3,3);
temp_power_zscore=zscore(temp_power).^2;
hold on;
for nCond=1:2
    for nCond2=1:2
        if nCond2==1
        temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Memory(~isnan(All_Conds))<median(All_Memory(~isnan(All_Conds))),:);
        else
        temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Memory(~isnan(All_Conds))>median(All_Memory(~isnan(All_Conds))),:);
        end
        if nCond2==1
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(myfreqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 30])
% ylim([-1 4])

hold on;
scatter(myfreqs(find(realpos_quad{1}.clusters)),-1+0.2*ones(1,length(find(realpos_quad{1}.clusters))),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(myfreqs(find(realpos_quad{2}.clusters)),-1.5+0.2*ones(1,length(find(realpos_quad{2}.clusters))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(myfreqs(find(realpos_quad{3}.clusters)),-2+0.2*ones(1,length(find(realpos_quad{3}.clusters))),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');


scatter(myfreqs(find(realneg_quad{1}.clusters)),-1.2+0.2*ones(1,length(find(realneg_quad{1}.clusters))),'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(myfreqs(find(realneg_quad{2}.clusters)),-1.7+0.2*ones(1,length(find(realneg_quad{2}.clusters))),'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(myfreqs(find(realneg_quad{3}.clusters)),-2+2.2*ones(1,length(find(realneg_quad{3}.clusters))),'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r');
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

% Freqs=[3.2 4.4; 9 9.8];
% Freqs=[1.4 3.4; 9 9.8];
Freqs=[7.8 10.6; 16.4 18.2]; % A modifier

figure;
hold on;
% cmap=cbrewer('seq','Blues',16);
% ylims=[-0.4 1.2];
% for nFreq=1:size(Freqs,1)
%     patch([Freqs(nFreq,1) Freqs(nFreq,2) Freqs(nFreq,2) Freqs(nFreq,1) Freqs(nFreq,1)],[ylims(1) ylims(1) ylims(2) ylims(2) ylims(1)],...
%         cmap(nFreq*4,:),'EdgeColor',cmap(nFreq*4,:),'FaceAlpha',0.5);
% end

set(gcf,'Position',[1     1   460   804]);
hp=[];
for nCond=1:2
    for nCond2=1:2
        temp_toplot=squeeze(Pow_allE(ismember(All_Conds,myConds{nCond}) & All_Memory==(nCond2-1),3,:));
        if nCond2==1
            hp(end+1)=plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            hp(end+1)=plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 30])
ylabel('Power (dB)')
%

scatter(myfreqs(find(realpos_quad2{2}.clusters2)),-0.2*ones(1,length(find(realpos_quad2{2}.clusters2))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[0 1 0]*0.8,'SizeData',64);
scatter(myfreqs(find(realneg_quad2{2}.clusters2)),-0.2*ones(1,length(find(realneg_quad2{2}.clusters2))),'Marker','s','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]*0.8,'SizeData',64);

%%
hl=legend(hp,{'WK | NO', 'WK | YS','SL | NO','SL | YE'},'Position',[0.7174    0.6405    0.2293    0.1063]);

% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_TimePlot.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_TimePlot.png'],'-r 300')

% Freqs=[1 4; 4 8; 8 12 ; 15 20];
TitleFreqs={'Delta','Alpha'};
figure;
set(gcf,'Position',[462     1   348   804]);
for nFreq=1
    subplot(size(Freqs,1),1,1);
    format_fig;
    for nCond=1:length(myConds)
        for nCond2=1:2
            temp_toplot=squeeze(Pow_allE(:,3,:));
            %         temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
            temp_toplot=mean(temp_toplot(:,freqs>=Freqs(nFreq,1) & freqs<=Freqs(nFreq,2)),2);
            temp_toplot=(nanzscore(temp_toplot));
            temp_toplot=temp_toplot(ismember(All_Conds,myConds{nCond}) & All_Memory==(nCond2-1),:);
            %%%%% Plot
            hold on;
            Pos=2*nCond2+0.2*(2*nCond-3); data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
            if nCond2==1
                colorBar=[ColorsGroup(nCond,:) ;ColorsGroup(nCond,:)];
            else
                colorBar=[ColorsGroup(nCond,:) ; 0 0 0];
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
    set(gca,'XTick',[2 4],'XTickLabel',{'NO','YES'})
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
            temp_toplot=temp_toplot(ismember(All_Conds,myConds{nCond}) & All_Memory==(nCond2-1),:);
            %%%%% Plot
            hold on;
            Pos=2*nCond2+0.2*(2*nCond-3); data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
            if nCond2==1
                colorBar=[ColorsGroup(nCond,:) ;ColorsGroup(nCond,:)];
            else
                colorBar=[ColorsGroup(nCond,:) ; 0 0 0];
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
    set(gca,'XTick',[2 4],'XTickLabel',{'NO','YES'})
    ylabel(sprintf('z(Power %s)^2',TitleFreqs{nFreq}))
    xlabel('Insight')
end
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_DotPlot.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_DotPlot.png'],'-r 300')

%%
figure;
stepbin=33;
clear bin_*
titlePlots={'Alpha','Sigma'};
set(gcf,'Position',[462     1   348   804]);
for nplot=1:3
    subplot(3,1,nplot);
    hold on;
    format_fig;
    if nplot==1
        tempX=nanzscore(Data.PowDelta);
    elseif nplot==2
        tempX=nanzscore(Data.PowAlpha);
    elseif nplot ==3
        tempX=nanzscore(Data.PowBeta);
        
    end
    %     tempY=Data.RT_fromAha-Data.RT_beforeAha;
    tempY=Data.Corrprenotpost;
    
    bins=prctile(tempX,0:stepbin:100);
    bin_values=[];
    bin_values{1,1}=(tempX(tempX<bins(2)));
    bin_values{1,2}=(tempY(tempX<bins(2)));
    for k=2:length(bins)-2
        bin_values{k,1}=(tempX(tempX>=bins(k) & tempX<=bins(k+1)));
        bin_values{k,2}=(tempY(tempX>=bins(k) & tempX<=bins(k+1)));
    end
    bin_values{length(bins)-1,1}=(tempX(tempX>bins(length(bins)-1)));
    bin_values{length(bins)-1,2}=(tempY(tempX>bins(length(bins)-1)));
    
    bin_values2=[];
    bin_values_sem=[];
    for kbin=1:size(bin_values,1)
        Pos=kbin; data=bin_values{kbin,2}; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
        colorBar=[.5 .5 .5; 0 0 0];
        
        
        xspread=(rand(1,length(data))-0.5)*widthBar/5+Pos;
        yspread=(rand(1,length(data))-0.5)*widthBar/10+data';
        scatter(xspread,yspread,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
        
        bin_values2(kbin,1)=mean(data);
        bin_values_sem(kbin,1)=sem(data);
    end
    errorbar(1:length(bins)-1,bin_values2(:,1),-bin_values_sem(:,1),+bin_values_sem(:,1),'Color','k','LineWidth',3);
    scatter(1:length(bins)-1,bin_values2(:,1),'Marker','o','SizeData',244,'MarkerFaceColor',[1 1 1]*0.7,'MarkerEdgeAlpha',.7,'MarkerEdgeColor','k','LineWidth',3);
    
    %     ylim([0 0.8])
    xlim([0.5 length(bins)-1+0.5])
    %     ylim([-0.08 1.08])
    set(gca,'XTick',1:length(bins),'XColor','k','YColor','k'); %,'XTickLabel',{'low','med','high'});
    %     title(titlePlots{nplot})
    ylabel('Insight','Color','k')
    xlabel({'Power Bin',titlePlots{nplot}},'Color','k')
end

% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_BinPlot.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_InsightEffect_BinPlot.png'],'-r 300')


%%
anovan(nanzscore(Data.PowAlpha),[Data.SleepNap Data.Corrprenotpost],'varnames',{'Sleep','Insight'},'model','full');
anovan(nanzscore(Data.PowSigma).^2,[Data.SleepNap Data.Corrprenotpost],'varnames',{'Sleep','Insight'},'model','full');

