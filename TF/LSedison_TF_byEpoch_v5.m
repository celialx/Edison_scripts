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

ColorsGroup=[178,171,210;
    253,184,99;
    230,97,1]/256;
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
all_AlphaPeaks=[];
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
    nc=nc+1;
    
    threshArt=150;
    fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(maxAbs(:,2)>threshArt),100*mean(maxAbs(:,2)>threshArt),threshArt)
    pption_rejTr(nc,:)=[sum(maxAbs(:,2)>threshArt) size(maxAbs,1) 100*mean(maxAbs(:,2)>threshArt)];
    for nCh=1:length(data.label)
        for nTr=1:length(data.trial)
            
            w_window=6*data.fsample;
            w_overlap=w_window/2;
            df=0.2;
            freqV=1:0.2:30;
            signal=data.trial{nTr}(nCh,:);
            [pow,faxis] = pwelch(signal,w_window,w_overlap,freqV,data.fsample,'psd');
            
            %             [faxis,pow]=get_PowerSpec(data.trial{nTr}(nCh,:), data.fsample, 0 ,0);
            TFRhann.powspctrm(nTr,nCh,:,1)=pow;
        end
    end
    
    AASM_Score=T.Stade;
    Beg_Task=(T.Epoch(find(T.Start(:,1))));
    End_Task=T.Epoch(find(T.End(:,1)));
    Drops=find(T.Ball(:,1));
    if ~isempty(Drops)
%         AASM_Score(Drops(1):end)=NaN;
    end
    AASM_Score=AASM_Score(T.Epoch>=Beg_Task & T.Epoch<=End_Task);
    AASM_Score_Epochs=(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score_Epochs2=unique(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score2=nan(1,length(data.trial));
    for j=1:length(data.trial)
        if sum(isnan(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)))))==length(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j))))
            AASM_Score2(j)=nan;
        else
            temp=unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)));
            AASM_Score2(j)=temp(1);
        end
    end
    AASM_Score2(find(maxAbs(:,3)>threshArt))=NaN;
    
    freqs=faxis;
    
    Stages=0:2;
    fooof_results_perStage=[];
    for nSta=1:4
        if nSta<4
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
            if sum(isnan(temp))==length(squeeze(temp)')
                fooof_results_perStage{nSta}.background_params=nan(1,2);
                
                FOOOF_Pow_perStage(nc,nSta,:)= nan(1,141);
                FOOOF2_Pow_perStage(nc,nSta,:)= nan(1,141);

            else
                fooof_results_perStage{nSta} = fooof(freqs,temp', f_range, settings,1);
                
                FOOOF_Pow_perStage(nc,nSta,:)= fooof_results_perStage{nSta}.fooofed_spectrum;
                FOOOF2_Pow_perStage(nc,nSta,:)= fooof_results_perStage{nSta}.fooofed_spectrum- fooof_results_perStage{nSta}.ap_fit;
            end
            Pow_perStage(nc,nSta,:)=squeeze(nanmean(mean(log(TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
        else
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(AASM_Score2~=0,2,:,:)),1),4));
            if sum(isnan(temp))==length(squeeze(temp)')
                fooof_results_perStage{nSta}.background_params=nan(1,2);
                
                FOOOF_Pow_perStage(nc,nSta,:)= nan(1,141);
                FOOOF2_Pow_perStage(nc,nSta,:)= nan(1,141);
            else
                fooof_results_perStage{nSta} = fooof(freqs,temp', f_range, settings,1);
            
            
                FOOOF_Pow_perStage(nc,nSta,:)= fooof_results_perStage{nSta}.fooofed_spectrum;
                FOOOF2_Pow_perStage(nc,nSta,:)= fooof_results_perStage{nSta}.fooofed_spectrum- fooof_results_perStage{nSta}.ap_fit;
           end
            Pow_perStage(nc,nSta,:)=squeeze(nanmean(mean(log(TFRhann.powspctrm(AASM_Score2~=0,2,:,:)),1),4));
        end
    end
    Pow_all(nc,:)=squeeze(nanmean(nanmean(log(TFRhann.powspctrm(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')),3,:,:)),1),4));
    for nCh=1:3
        Pow_allE(nc,nCh,:)=squeeze(nanmean(nanmean(log(TFRhann.powspctrm(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')),nCh,:,:)),1),4));
        fooof_results = fooof(freqs,squeeze(nanmean(nanmean((TFRhann.powspctrm(find(maxAbs(:,nCh)<=threshArt & ~isnan(AASM_Score2')),nCh,:,:)),1),4))', f_range, settings,1);
        FOOOF_Pow_allE(nc,nCh,:)=fooof_results.fooofed_spectrum;
        FOOOF2_Pow_allE(nc,nCh,:)=fooof_results.fooofed_spectrum-fooof_results.ap_fit;
    end
    fooof_results = fooof(freqs,squeeze(nanmean(nanmean((TFRhann.powspctrm(find(maxAbs(:,3)<=threshArt & ~isnan(AASM_Score2')),3,:,:)),1),4))', f_range, settings,1);

    all_Pow=[all_Pow ; [nF mean(Pow_all(nc,freqs>2 & freqs<4)) mean(Pow_all(nc,freqs>4 & freqs<8)) mean(Pow_all(nc,freqs>8 & freqs<11)) mean(Pow_all(nc,freqs>1 & freqs<40)) fooof_results.background_params length(AASM_Score2) sum(AASM_Score2==0) sum(AASM_Score2==1) sum(AASM_Score2==2) ...
        mean((Pow_perStage(nc,1,freqs>2 & freqs<4))) mean((Pow_perStage(nc,1,freqs>4 & freqs<8))) mean((Pow_perStage(nc,1,freqs>8 & freqs<11))) mean((Pow_perStage(nc,1,freqs>1 & freqs<40))) fooof_results_perStage{1}.background_params...
        mean((Pow_perStage(nc,2,freqs>2 & freqs<4))) mean((Pow_perStage(nc,2,freqs>4 & freqs<8))) mean((Pow_perStage(nc,2,freqs>8 & freqs<11))) mean((Pow_perStage(nc,2,freqs>1 & freqs<40))) fooof_results_perStage{2}.background_params...
        mean((Pow_perStage(nc,3,freqs>2 & freqs<4))) mean((Pow_perStage(nc,3,freqs>4 & freqs<8))) mean((Pow_perStage(nc,3,freqs>8 & freqs<11))) mean((Pow_perStage(nc,3,freqs>1 & freqs<40))) fooof_results_perStage{3}.background_params...
        mean((Pow_perStage(nc,4,freqs>2 & freqs<4))) mean((Pow_perStage(nc,4,freqs>4 & freqs<8))) mean((Pow_perStage(nc,4,freqs>8 & freqs<11))) mean((Pow_perStage(nc,4,freqs>1 & freqs<40))) fooof_results_perStage{4}.background_params ]];
    all_Subds=[all_Subds ; {SubdID}];
    
    alpha_peaks=fooof_results.peak_params(fooof_results.peak_params(:,1)>7 & fooof_results.peak_params(:,1)<12,:);
    if ~isempty(alpha_peaks)
        alpha_peaks=alpha_peaks(alpha_peaks(:,2)==max(alpha_peaks(:,2)),:);
         all_AlphaPeaks=[all_AlphaPeaks ; [nF sum(AASM_Score==0) sum(AASM_Score==1) sum(AASM_Score==2) alpha_peaks]];
    else
         all_AlphaPeaks=[all_AlphaPeaks ; [nF sum(AASM_Score==0) sum(AASM_Score==1) sum(AASM_Score==2) nan(1,3)]];
    end
    %     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
end
freqs2=fooof_results.freqs;

%%
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
    
    data_clean.PowDeltaTheta(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>1.2 & freqs2<7.6),3));
    data_clean.PowDelta(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>1 & freqs2<4),3));
    data_clean.PowTheta(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>4 & freqs2<8),3));
    data_clean.PowAlpha(this_line)=squeeze(mean(FOOOF_Pow_allE(nS,3,freqs2>8 & freqs2<12),3));
    
    data_clean.PowDelta_dev(this_line)=squeeze(mean(FOOOF2_Pow_allE(nS,3,freqs2>1 & freqs2<4),3));
    data_clean.PowTheta_dev(this_line)=squeeze(mean(FOOOF2_Pow_allE(nS,3,freqs2>4 & freqs2<8),3));
    data_clean.PowAlpha_dev(this_line)=squeeze(mean(FOOOF2_Pow_allE(nS,3,freqs2>8 & freqs2<12),3));
    
    
    data_clean.PowDeltaN(this_line)=all_Pow(nS,2)-all_Pow(nS,5);
    data_clean.PowThetaN(this_line)=all_Pow(nS,3)-all_Pow(nS,5);
    data_clean.PowAlphaN(this_line)=all_Pow(nS,4)-all_Pow(nS,5);
    
    data_clean.PowDeltaN_W(this_line)=all_Pow(nS,12)-all_Pow(nS,15);
    data_clean.PowThetaN_W(this_line)=all_Pow(nS,13)-all_Pow(nS,15);
    data_clean.PowAlphaN_W(this_line)=all_Pow(nS,14)-all_Pow(nS,15);
    
    data_clean.PowDeltaN_S(this_line)=all_Pow(nS,30)-all_Pow(nS,33);
    data_clean.PowThetaN_S(this_line)=all_Pow(nS,31)-all_Pow(nS,33);
    data_clean.PowAlphaN_S(this_line)=all_Pow(nS,32)-all_Pow(nS,33);
    
    data_clean.PowThetaAlpha(this_line)=mean(Pow_all(nS,freqs>2 & freqs<7))./mean(Pow_all(nS,freqs>7 & freqs<12));
    
    data_clean.PowThetaDelta(this_line)=mean(Pow_all(nS,freqs>5 & freqs<8))./mean(Pow_all(nS,freqs>2 & freqs<4));
    
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

% writetable(data_clean,'../Edison_Tables/Clean_Data_SS_Pow.txt');
%%
figure; hold on;
for nCond=1:3
    temp_toplot=squeeze(FOOOF_Pow_perStage(:,nCond,:));
    %     temp_toplot=(temp_toplot)-repmat(nanmean(temp_toplot,2),1,length(freqs));
    plot(freqs2,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
end
xlabel('Frequency (Hz)')
legend({'WK','N1','N2'})
format_fig;
xlim([2 20])
% ylim([-1 3])
%%
myConds={[1],[2 3]};
freqs=freqs;
figure; hold on;
for nCond=1:3
    temp_toplot=squeeze(FOOOF_Pow_allE(All_Conds==nCond,3,:));
    %     temp_toplot=(temp_toplot)-repmat(nanmean(temp_toplot,2),1,length(freqs));
    plot(freqs2,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
end
xlabel('Frequency (Hz)')
legend({'WK','N1','N2'})
format_fig;
xlim([2 20])
% ylim([-1 3])



%%
Freqs=[1 4; 4 8; 8 12 ; 15 20];
TitleFreqs={'detla','theta','alpha','sigma'};

for nFreq=1:size(Freqs,1)
    figure;
    set(gcf,'Position',[287   325   291   419]);
    format_fig;
    for nCond=1:length(myConds)
        for nCond2=1:2
            temp_toplot=squeeze(FOOOF_Pow_allE(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),2,:));
            %         temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
            temp_toplot=mean(temp_toplot(:,freqs>Freqs(nFreq,1) & freqs<Freqs(nFreq,2)),2);
            
            %%%%% Plot
            hold on;
            Pos=2*nCond+0.2*(2*nCond2-3); data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
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
            
            xspread=(rand(1,length(data))-0.5)*widthBar/4+2*nCond+.65*(2*nCond2-3);
            scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
            
            
            %         [hdot pV]=simpleDotPlot(nCond+0.2*(2*nCond2-3),temp_toplot,400,ColorsGroup(nCond,:),1,ColorsGroup(nCond,:),'o',[],3,1,1);
            %         if nCond2==2
            %             set(hdot.ind,'MarkerEdgeColor','k');
            %         end
        end
        %     hold on;
        %     scatter(hdot.ind.XData(temp_insight==1),hdot.ind.YData(temp_insight==1),...
        %         'Marker','o','MarkerFaceColor',ColorsGroup(nCond,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'SizeData',100,'LineWidth',2);
    end
    % xlim([0.2 3.8])
    % set(gca,'XTick',1:3,'XTickLabel',{'WK','N1','N2'})
    xlim([.5 5.5])
    set(gca,'XTick',[2 4],'XTickLabel',{'WK','SL'})
    title(TitleFreqs{nFreq})
    ylabel('Power')
end

%% STATS
All_Conds2=All_Conds~=1;
clear Effect*
for nF=1:size(FOOOF_Pow_allE,3)
    [p,t,stats,terms] = anovan(FOOOF_Pow_allE((~isnan(All_Conds)),3,nF),[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],'model','full','display','off');
    EffectSleep(nF,1)=t{2,6};
    EffectInsigth(nF,1)=t{3,6};
    EffectInteraction(nF,1)=t{4,6};
    
    EffectSleep(nF,2)=t{2,7};
    EffectInsigth(nF,2)=t{3,7};
    EffectInteraction(nF,2)=t{4,7};
end


figure;
plot(freqs2,EffectSleep(:,1),'k');
hold on;
plot(freqs2,EffectInsigth(:,1),'r');
plot(freqs2,EffectInteraction(:,1),'b');

%%

% [realpos realneg]=get_cluster_permutation_aov(squeeze(Pow_allE((~isnan(All_Conds)),3,:)),[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],...
%     0.05,0.1,100,freqs);

temp_power=squeeze(Pow_allE((~isnan(All_Conds)),3,:));
[realpos_lin realneg_lin]=get_cluster_permutation_aov(zscore(temp_power),[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],...
    0.05,0.1,1000,freqs);
[realpos_quad realneg_quad]=get_cluster_permutation_aov(zscore(temp_power).^2,[All_Conds2(~isnan(All_Conds)) All_Insight(~isnan(All_Conds))],...
    0.05,0.1,1000,freqs);
%%
figure; 
subplot(1,3,1);
hold on;
for nCond=1:2
    for nCond2=1:2
        temp_toplot=squeeze(Pow_allE(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),3,:));
        if nCond2==1
            plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
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
            plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 20])
% ylim([-1 4])

hold on;
scatter(freqs(find(realpos_lin{1}.clusters)),-1+0.2*ones(1,length(find(realpos_lin{1}.clusters))),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(freqs(find(realpos_lin{2}.clusters)),-1.5+0.2*ones(1,length(find(realpos_lin{2}.clusters))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(freqs(find(realpos_lin{3}.clusters)),-2+0.2*ones(1,length(find(realpos_lin{3}.clusters))),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');


subplot(1,3,3);
temp_power_zscore=zscore(temp_power).^2;
hold on;
for nCond=1:2
    for nCond2=1:2
        temp_toplot=temp_power_zscore(ismember(All_Conds(~isnan(All_Conds)),myConds{nCond}) & All_Insight(~isnan(All_Conds))==(nCond2-1),:);
        if nCond2==1
            plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
            plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
format_fig;
xlim([1 20])
% ylim([-1 4])

hold on;
scatter(freqs(find(realpos_quad{1}.clusters)),-1+0.2*ones(1,length(find(realpos_quad{1}.clusters))),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b');
scatter(freqs(find(realpos_quad{2}.clusters)),-1.5+0.2*ones(1,length(find(realpos_quad{2}.clusters))),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(freqs(find(realpos_quad{3}.clusters)),-2+0.2*ones(1,length(find(realpos_quad{3}.clusters))),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');



%%
mdl = fitglm(data_clean,'InsightPost~1+SleepEdison + PowDeltaTheta','Distribution','binomial');
writetable(data_clean,'Edison_Table_Pow.csv');

figure; set(gcf,'Position',[440   347   388   451]);
bar(1,mdl.Coefficients.Estimate(1),'FaceColor',[0.6 0.6 .6],'EdgeColor',[0 0 0],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(1,mdl.Coefficients.Estimate(1),mdl.Coefficients.SE(1),'Color',[0 0 0],'LineWidth',3);


bar(2,mdl.Coefficients.Estimate(2),'FaceColor',[0 0 .6],'EdgeColor',[0 0 1],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(2,mdl.Coefficients.Estimate(2),mdl.Coefficients.SE(2),'Color',[0 0 1],'LineWidth',3);

bar(3,mdl.Coefficients.Estimate(3),'FaceColor',[0.6 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(3,mdl.Coefficients.Estimate(3),mdl.Coefficients.SE(3),'Color',[1 0 0],'LineWidth',3);


xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'Int','Sleep','\delta~\theta'});
format_fig;
ylabel('Coefficient')

%%
mdl = fitglm(data_clean,'InsightPost~1+SleepGroup + PowDeltaTheta','Distribution','binomial');

figure; set(gcf,'Position',[440   347   471   451]);
bar(1,mdl.Coefficients.Estimate(1),'FaceColor',[0.6 0.6 .6],'EdgeColor',[0 0 0],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(1,mdl.Coefficients.Estimate(1),mdl.Coefficients.SE(1),'Color',[0 0 0],'LineWidth',3);


bar(2,mdl.Coefficients.Estimate(2),'FaceColor',[0 0 .6],'EdgeColor',[0 0 1],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(2,mdl.Coefficients.Estimate(2),mdl.Coefficients.SE(2),'Color',[0 0 1],'LineWidth',3);

bar(3,mdl.Coefficients.Estimate(3),'FaceColor',[1 1 1],'EdgeColor',[0 0 1],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(3,mdl.Coefficients.Estimate(3),mdl.Coefficients.SE(3),'Color',[0 0 1],'LineWidth',3);

bar(4,mdl.Coefficients.Estimate(4),'FaceColor',[0.6 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.5,'LineWidth',3);
hold on;
errorbar(4,mdl.Coefficients.Estimate(4),mdl.Coefficients.SE(4),'Color',[1 0 0],'LineWidth',3);


xlim([0.2 4.8])
set(gca,'XTick',1:4,'XTickLabel',{'Int','N1','N2','\delta~\theta'});
format_fig;
ylabel('Coefficient')


%%
stepbin=33;
clear bin_*
titlePlots={'\delta','\theta','\alpha','\alpha peak'};
figure; set(gcf,'Position',[199   421   840   377]);
for nplot=1:3
    subplot(1,3,nplot);
    hold on;
    format_fig;
    if nplot==1
        tempX=nanzscore(data_clean.PowDelta);
    elseif nplot==2
        tempX=nanzscore(data_clean.PowTheta);
    elseif nplot==3
        tempX=nanzscore(data_clean.PowAlpha);
%         elseif nplot==4
%         tempX=nanzscore(data_clean.PowAlpha); %
%         tempX=(minmax(nanzscore(data_clean.PowAlpha))+1)./(minmax(nanzscore(data_clean.PowTheta))+1);
    end
%     tempY=data_clean.RT_fromAha-data_clean.RT_beforeAha;
    tempY=data_clean.InsightPost;
    
    bins=prctile(tempX,0:stepbin:100);
    bin_values(1,1)=mean(tempX(tempX<bins(2)));
    bin_values(1,2)=mean(tempY(tempX<bins(2)));
    for k=2:length(bins)-2
    bin_values(k,1)=mean(tempX(tempX>=bins(k) & tempX<=bins(k+1)));
    bin_values(k,2)=mean(tempY(tempX>=bins(k) & tempX<=bins(k+1)));
    end
    bin_values(length(bins)-1,1)=mean(tempX(tempX>bins(length(bins)-1)));
    bin_values(length(bins)-1,2)=mean(tempY(tempX>bins(length(bins)-1)));
    
    bin_values_sem(1,1)=sem(tempX(tempX<bins(2)));
    bin_values_sem(1,2)=sem(tempY(tempX<bins(2)));
    for k=2:length(bins)-2
    bin_values_sem(k,1)=sem(tempX(tempX>=bins(k) & tempX<=bins(k+1)));
    bin_values_sem(k,2)=sem(tempY(tempX>=bins(k) & tempX<=bins(k+1)));
    end
    bin_values_sem(length(bins)-1,1)=sem(tempX(tempX>bins(length(bins)-1)));
    bin_values_sem(length(bins)-1,2)=sem(tempY(tempX>bins(length(bins)-1)));
    
    errorbar(1:length(bins)-1,bin_values(:,2),-bin_values_sem(:,2)/2,+bin_values_sem(:,2)/2,'Color','k','LineWidth',3);
    scatter(1:length(bins)-1,bin_values(:,2),'Marker','o','SizeData',288,'MarkerFaceColor',[1 1 1]*0.7,'MarkerEdgeAlpha',.7,'MarkerEdgeColor','k','LineWidth',3);
%     ylim([0 0.8])
    xlim([0.5 length(bins)-1+0.5])
    set(gca,'XTick',1:length(bins)); %,'XTickLabel',{'low','med','high'});
    title(titlePlots{nplot})
    ylabel('RT diff')
    xlabel('Power Bin')
end

%%
%%
titlePlots={'offset','slope'};
figure; set(gcf,'Position',[199   416   801   382]);
for nplot=1:2
    subplot(1,2,nplot);
    hold on;
    format_fig;
    if nplot==1
        tempX=nanzscore(data_clean.PowBG);
    elseif nplot==2
        tempX=nanzscore(data_clean.PowSlope);
    end
    tempY=data_clean.InsightPost;
    
    bins=prctile(tempX,0:stepbin:100);
    bin_values(1,1)=mean(tempX(tempX<bins(2)));
    bin_values(1,2)=mean(tempY(tempX<bins(2)));
    for k=2:length(bins)-2
    bin_values(k,1)=mean(tempX(tempX>=bins(k) & tempX<=bins(k+1)));
    bin_values(k,2)=mean(tempY(tempX>=bins(k) & tempX<=bins(k+1)));
    end
    bin_values(length(bins)-1,1)=mean(tempX(tempX>bins(length(bins)-1)));
    bin_values(length(bins)-1,2)=mean(tempY(tempX>bins(length(bins)-1)));
    
    bin_values_sem(1,1)=sem(tempX(tempX<bins(2)));
    bin_values_sem(1,2)=sem(tempY(tempX<bins(2)));
    for k=2:length(bins)-2
    bin_values_sem(k,1)=sem(tempX(tempX>=bins(k) & tempX<=bins(k+1)));
    bin_values_sem(k,2)=sem(tempY(tempX>=bins(k) & tempX<=bins(k+1)));
    end
    bin_values_sem(length(bins)-1,1)=sem(tempX(tempX>bins(length(bins)-1)));
    bin_values_sem(length(bins)-1,2)=sem(tempY(tempX>bins(length(bins)-1)));
    
    errorbar(1:length(bins)-1,bin_values(:,2),-bin_values_sem(:,2)/2,+bin_values_sem(:,2)/2,'Color','k','LineWidth',3);
    scatter(1:length(bins)-1,bin_values(:,2),'Marker','o','SizeData',288,'MarkerFaceColor',[1 1 1]*0.7,'MarkerEdgeAlpha',.7,'MarkerEdgeColor','k','LineWidth',3);
    ylim([0.1 0.75])
    xlim([0.5 length(bins)-1+0.5])
    set(gca,'XTick',1:length(bins)); %,'XTickLabel',{'low','med','high'});
    title(titlePlots{nplot})
    ylabel('Insight Prob')
    xlabel('Bin')
end

%%
bins=prctile(data_clean.PowDelta,0:33:100);  data_clean.PowDeltaBin=data_clean.PowDelta;
data_clean.PowDeltaBin(data_clean.PowDelta<bins(2))=1;
data_clean.PowDeltaBin(data_clean.PowDelta>=bins(2) & data_clean.PowDelta<=bins(3))=2;
data_clean.PowDeltaBin(data_clean.PowDelta>bins(3))=3;
bins=prctile(data_clean.PowTheta,0:33:100); data_clean.PowThetaBin=data_clean.PowTheta;
data_clean.PowThetaBin(data_clean.PowTheta<bins(2))=1;
data_clean.PowThetaBin(data_clean.PowTheta>=bins(2) & data_clean.PowTheta<=bins(3))=2;
data_clean.PowThetaBin(data_clean.PowTheta>bins(3))=3;
bins=prctile(data_clean.PowAlpha,0:33:100); data_clean.PowAlphaBin=data_clean.PowAlpha;
data_clean.PowAlphaBin(data_clean.PowAlpha<bins(2))=1;
data_clean.PowAlphaBin(data_clean.PowAlpha>=bins(2) & data_clean.PowAlpha<=bins(3))=2;
data_clean.PowAlphaBin(data_clean.PowAlpha>bins(3))=3;

bins=prctile(data_clean.PowDelta_dev,0:33:100);  data_clean.PowDeltaBin_dev=data_clean.PowDelta_dev;
data_clean.PowDeltaBin_dev(data_clean.PowDelta_dev<bins(2))=1;
data_clean.PowDeltaBin_dev(data_clean.PowDelta_dev>=bins(2) & data_clean.PowDelta_dev<=bins(3))=2;
data_clean.PowDeltaBin_dev(data_clean.PowDelta_dev>bins(3))=3;
bins=prctile(data_clean.PowTheta_dev,0:33:100); data_clean.PowThetaBin_dev=data_clean.PowTheta_dev;
data_clean.PowThetaBin_dev(data_clean.PowTheta_dev<bins(2))=1;
data_clean.PowThetaBin_dev(data_clean.PowTheta_dev>=bins(2) & data_clean.PowTheta_dev<=bins(3))=2;
data_clean.PowThetaBin_dev(data_clean.PowTheta_dev>bins(3))=3;
bins=prctile(data_clean.PowAlpha_dev,0:33:100); data_clean.PowAlphaBin_dev=data_clean.PowAlpha_dev;
data_clean.PowAlphaBin_dev(data_clean.PowAlpha_dev<bins(2))=1;
data_clean.PowAlphaBin_dev(data_clean.PowAlpha_dev>=bins(2) & data_clean.PowAlpha_dev<=bins(3))=2;
data_clean.PowAlphaBin_dev(data_clean.PowAlpha_dev>bins(3))=3;

bins=prctile(data_clean.PowBG,0:33:100);  data_clean.PowBGBin=data_clean.PowBG;
data_clean.PowBGBin(data_clean.PowBG<bins(2))=1;
data_clean.PowBGBin(data_clean.PowBG>=bins(2) & data_clean.PowBG<=bins(3))=2;
data_clean.PowBGBin(data_clean.PowBG>bins(3))=3;
bins=prctile(data_clean.PowSlope,0:33:100); data_clean.PowSlopeBin=data_clean.PowSlope;
data_clean.PowSlopeBin(data_clean.PowSlope<bins(2))=1;
data_clean.PowSlopeBin(data_clean.PowSlope>=bins(2) & data_clean.PowSlope<=bins(3))=2;
data_clean.PowSlopeBin(data_clean.PowSlope>bins(3))=3;

writetable(data_clean,'../Edison_Tables/Clean_Data_SS_Pow_beforeDrop.txt');

%%
tempY=data_clean.RT_fromAha-data_clean.RT_beforeAha;
tempX=data_clean.PowDelta;
[p0,S0,mu0] = polyfit(tempX(~isnan(tempX)),tempY(~isnan(tempX)),0);
[p1,S1,mu1] = polyfit(tempX(~isnan(tempX)),tempY(~isnan(tempX)),1);
[p2,S2,mu2] = polyfit(tempX(~isnan(tempX)),tempY(~isnan(tempX)),2);


%%
% figure; set(gcf,'Position',[ 440   303   535   495]);
% format_fig;
% cmap=cbrewer('seq','Blues',4);
% simpleBarPlot(1,100*mean(data_clean.InsightPost(data_clean.SleepGroup=='0')),cmap(2,:),0.85,'k');
% simpleBarPlot(2,100*mean(data_clean.InsightPost(data_clean.SleepGroup=='1')),cmap(3,:),0.85,'k');
% simpleBarPlot(3,100*mean(data_clean.InsightPost(data_clean.SleepGroup=='2')),cmap(4,:),0.85,'k');
% xlim([0.2 3.8])
% set(gca,'XTick',1:3,'XTickLabel',{'WAKE','N1','N2'},'FontSize',22);
% ylabel('Insight (%)')
% xlabel('Sleep Group')