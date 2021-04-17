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
    fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(maxAbs(:,2)>threshArt),100*mean(maxAbs(:,2)>threshArt),threshArt)
    
    cfg              = [];
    cfg.trials       = find(maxAbs(:,2)<=threshArt);
    cfg.output       = 'pow';
    cfg.channel      = 'all';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
    cfg.toi          = 'all';                         % time
    cfg.keeptrials  = 'yes';
    TFRhann = ft_freqanalysis(cfg, data);
    
    
    AASM_Score=T.Stade;
    Beg_Task=(T.Epoch(find(T.Start(:,1))));
    End_Task=T.Epoch(find(T.End(:,1)));
    AASM_Score=AASM_Score(T.Epoch>=Beg_Task & T.Epoch<=End_Task);
    AASM_Score_Epochs=(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score_Epochs2=unique(T.Epoch(T.Epoch>=Beg_Task & T.Epoch<=End_Task));
    AASM_Score2=nan(1,length(data.trial));
    for j=1:length(data.trial)
        if sum(isnan(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)))))==length(unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j))))
            AASM_Score2(j)=nan;
        else
            AASM_Score2(j)=unique(AASM_Score(AASM_Score_Epochs==AASM_Score_Epochs2(j)));
        end
    end
    AASM_Score2=AASM_Score2(find(maxAbs(:,2)<=threshArt));
    freqs=TFRhann.freq;
    
    Stages=0:2;
    fooof_results_perStage=[];
    for nSta=1:4
        if nSta<4
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
            if sum(isnan(temp))==length(squeeze(temp)')
                fooof_results_perStage{nSta}.background_params=nan(1,2);
            else
                fooof_results_perStage{nSta} = fooof(freqs,temp', f_range, settings,1);
            end
            Pow_perStage(nc,nSta,:)=squeeze(nanmean(mean(log(TFRhann.powspctrm(AASM_Score2==Stages(nSta),2,:,:)),1),4));
        else
            temp=squeeze(nanmean(mean((TFRhann.powspctrm(AASM_Score2~=0,2,:,:)),1),4));
            if sum(isnan(temp))==length(squeeze(temp)')
                fooof_results_perStage{nSta}.background_params=nan(1,2);
            else
                fooof_results_perStage{nSta} = fooof(freqs,temp', f_range, settings,1);
            end
            Pow_perStage(nc,nSta,:)=squeeze(nanmean(mean(log(TFRhann.powspctrm(AASM_Score2~=0,2,:,:)),1),4));
        end
    end
    Pow_all(nc,:)=squeeze(nanmean(nanmean(log(TFRhann.powspctrm(:,2,:,:)),1),4));
    fooof_results = fooof(freqs,squeeze(nanmean(nanmean((TFRhann.powspctrm(:,2,:,:)),1),4))', f_range, settings,1);
    
    all_Pow=[all_Pow ; [nF mean(Pow_all(nc,freqs>2 & freqs<4)) mean(Pow_all(nc,freqs>4 & freqs<8)) mean(Pow_all(nc,freqs>8 & freqs<11)) mean(Pow_all(nc,freqs>1 & freqs<40)) fooof_results.background_params length(AASM_Score2) sum(AASM_Score2==0) sum(AASM_Score2==1) sum(AASM_Score2==2) ...
        mean((Pow_perStage(nc,1,freqs>2 & freqs<4))) mean((Pow_perStage(nc,1,freqs>4 & freqs<8))) mean((Pow_perStage(nc,1,freqs>8 & freqs<11))) mean((Pow_perStage(nc,1,freqs>1 & freqs<40))) fooof_results_perStage{1}.background_params...
        mean((Pow_perStage(nc,2,freqs>2 & freqs<4))) mean((Pow_perStage(nc,2,freqs>4 & freqs<8))) mean((Pow_perStage(nc,2,freqs>8 & freqs<11))) mean((Pow_perStage(nc,2,freqs>1 & freqs<40))) fooof_results_perStage{2}.background_params...
        mean((Pow_perStage(nc,3,freqs>2 & freqs<4))) mean((Pow_perStage(nc,3,freqs>4 & freqs<8))) mean((Pow_perStage(nc,3,freqs>8 & freqs<11))) mean((Pow_perStage(nc,3,freqs>1 & freqs<40))) fooof_results_perStage{3}.background_params...
        mean((Pow_perStage(nc,4,freqs>2 & freqs<4))) mean((Pow_perStage(nc,4,freqs>4 & freqs<8))) mean((Pow_perStage(nc,4,freqs>8 & freqs<11))) mean((Pow_perStage(nc,4,freqs>1 & freqs<40))) fooof_results_perStage{4}.background_params ]];
    all_Subds=[all_Subds ; {SubdID}];
    %     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
end

%%
data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');
data_clean.PowDelta=nan(size(data_clean,1),1);
data_clean.PowTheta=nan(size(data_clean,1),1);
data_clean.PowAlpha=nan(size(data_clean,1),1);

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

data_clean.PowSlope=nan(size(data_clean,1),1);
data_clean.PowBG=nan(size(data_clean,1),1);

data_clean.SleepGroup=nan(size(data_clean,1),1);
data_clean.SleepGroup(data_clean.SleepEdison==0)=0;
data_clean.SleepGroup(data_clean.SleepEdison==1  & data_clean.HypnoEdisonN2==0)=1;
data_clean.SleepGroup(data_clean.SleepEdison==1  & data_clean.HypnoEdisonN2~=0)=2;
data_clean.SleepGroup=categorical(data_clean.SleepGroup);
% data_clean.SleepGroup=reordercats(data_clean.SleepGroup,[2 1 3]);


All_Conds=[];
All_Insight=[];
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
    
    data_clean.PowDeltaN_S(this_line)=all_Pow(nS,30)-all_Pow(nS,33);
    data_clean.PowThetaN_S(this_line)=all_Pow(nS,31)-all_Pow(nS,33);
    data_clean.PowAlphaN_S(this_line)=all_Pow(nS,32)-all_Pow(nS,33);
    
    data_clean.PowThetaAlpha(this_line)=mean(Pow_all(nS,TFRhann.freq>2 & TFRhann.freq<7))./mean(Pow_all(nS,TFRhann.freq>7 & TFRhann.freq<12));

        data_clean.PowThetaDelta(this_line)=mean(Pow_all(nS,TFRhann.freq>5 & TFRhann.freq<8))./mean(Pow_all(nS,TFRhann.freq>2 & TFRhann.freq<4));

    data_clean.PowSlope(this_line)=all_Pow(nS,7);
    data_clean.PowBG(this_line)=all_Pow(nS,6);
    
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
freqs=TFRhann.freq;
times=TFRhann.time;
figure; hold on;
for nCond=1:3
    temp_toplot=squeeze(Pow_all(All_Conds==nCond,:));
    temp_toplot=(temp_toplot)-repmat(nanmean(temp_toplot,2),1,length(freqs));
    plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
end
xlabel('Frequency (Hz)')
legend({'WK','N1','N2'})
format_fig;
xlim([2 20])
ylim([-1 3])

%%
freqs=TFRhann.freq;
times=TFRhann.time;
figure; hold on;
for nCond=1:2
    for nCond2=1:2
        temp_toplot=squeeze(Pow_all(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),:));
%         temp_toplot=(temp_toplot)-repmat(nanmean(temp_toplot(:,freqs<20),2),1,length(freqs));
        if nCond2==1
        plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3,'LineStyle','--');
        else
        plot(freqs,nanmean(temp_toplot),'Color',ColorsGroup(nCond,:),'LineWidth',3);
        end
    end
end
xlabel('Frequency (Hz)')
legend({'WK | NO', 'WK | YS','SL | NO','SL | YE'})
format_fig;
xlim([1 20])
% ylim([-1 4])


%%
Freqs=[1 4; 4 7; 8 12];
TitleFreqs={'detla','theta','alpha'};
for nFreq=1:size(Freqs,1)
figure;
set(gcf,'Position',[287   325   291   419]);
format_fig;
myConds={[1],[2 3]};
for nCond=1:length(myConds)
    for nCond2=1:2
        temp_toplot=squeeze(Pow_all(ismember(All_Conds,myConds{nCond}) & All_Insight==(nCond2-1),:));
        temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
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
%%
figure;
set(gcf,'Position',[287   325   291   419]);
format_fig;
for nCond=1:3
    temp_toplot=squeeze(Pow_all(All_Conds==nCond,:));
%     temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
    temp_toplot=mean(temp_toplot(:,freqs>5 & freqs<8),2);
        temp_insight=squeeze(All_Insight(All_Conds==nCond));
[hdot pV]=simpleDotPlot(nCond,temp_toplot,400,ColorsGroup(nCond,:),2.5,ColorsGroup(nCond,:),'o',[],3,1,1);

    hold on;
    scatter(hdot.ind.XData(temp_insight==1),hdot.ind.YData(temp_insight==1),...
        'Marker','o','MarkerFaceColor',ColorsGroup(nCond,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'SizeData',100,'LineWidth',2);
end
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'WK','N1','N2'})
title('theta')

%%
figure;
set(gcf,'Position',[287   325   291   419]);
format_fig;
for nCond=1:3
    temp_toplot=squeeze(Pow_all(All_Conds==nCond,:));
    temp_toplot=(temp_toplot-repmat(nanmean(temp_toplot,2),1,length(freqs)));
    temp_toplot=mean(temp_toplot(:,freqs>8 & freqs<11),2);
        temp_insight=squeeze(All_Insight(All_Conds==nCond));
[hdot pV]=simpleDotPlot(nCond,temp_toplot,400,ColorsGroup(nCond,:),2.5,ColorsGroup(nCond,:),'o',[],3,1,1);

    hold on;
    scatter(hdot.ind.XData(temp_insight==1),hdot.ind.YData(temp_insight==1),...
        'Marker','o','MarkerFaceColor',ColorsGroup(nCond,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'SizeData',100,'LineWidth',2);
end
xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'WK','N1','N2'})
title('alpha')
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

%% Slept or Not - Delta Power
varname='PowDeltaN_S';
figure;
temp=data_clean.(varname)(data_clean.SleepEdison==0 & data_clean.InsightPost==0);
simpleBarPlot(1-0.2,temp,[1 1 1;1 0 0],0.35,'k');

temp=data_clean.(varname)(data_clean.SleepEdison==0 & data_clean.InsightPost==1);
simpleBarPlot(1+0.2,temp,[.7 .7 .7;1 0 0],0.35,'k');


temp=data_clean.(varname)(data_clean.SleepEdison==1  & data_clean.HypnoEdisonN2==0 & data_clean.InsightPost==0);
simpleBarPlot(2-0.2,temp,[1 1 1;0 0 1],0.35,'k');

temp=data_clean.(varname)(data_clean.SleepEdison==1  & data_clean.HypnoEdisonN2==0 & data_clean.InsightPost==1);
simpleBarPlot(2+0.2,temp,[.7 .7 .7;0 0 1],0.35,'k');

temp=data_clean.(varname)(data_clean.SleepEdison==1  & data_clean.HypnoEdisonN2~=0 & data_clean.InsightPost==0);
simpleBarPlot(3-0.2,temp,[1 1 1;0 0 1],0.35,'k');

temp=data_clean.(varname)(data_clean.SleepEdison==1  & data_clean.HypnoEdisonN2~=0 & data_clean.InsightPost==1);
simpleBarPlot(3+0.2,temp,[.7 .7 .7;0 0 1],0.35,'k');


xlim([0.2 3.8])
set(gca,'XTick',1:3,'XTickLabel',{'W','N1','N2'});
xlabel('Slept')
ylabel('Delta Power')
format_fig

%%

[B,DEV,STATS] = glmfit(data_clean.SleepEdison,data_clean.InsightPost,'binomial');
[B2,DEV2,STATS2] = glmfit([data_clean.SleepEdison data_clean.PowDelta],data_clean.InsightPost,'binomial');

%% Slept or Not - Delta Power
varname='PowTheta';
figure;
temp=data_clean.(varname)(data_clean.SleepGroup=='0');
simpleBarPlot(1,temp,[1 1 1;1 0 0],0.85,'k');

temp=data_clean.(varname)(data_clean.SleepGroup=='1');
simpleBarPlot(2,temp,[1 1 1;0 0 1],0.85,'k');

temp=data_clean.(varname)(data_clean.SleepGroup=='2');
simpleBarPlot(3,temp,[.7 .7 .7;0 0 1],0.85,'k');

xlim([0.2 3.8])
set(gca,'XTick',1:2,'XTickLabel',{'W','N1','N2'});
xlabel('Slept')
ylabel('Delta Power')
format_fig

%%
varname='PowDeltaN';
figure;
% bins=prctile(data_clean.(varname),0:20:100);
% thetaP=[];
% thetaN=[];
% thetaSEM=[];
% insightP=[];
% insightN=[];
% insightSEM=[];
% for i=1:length(bins)-1
%     thetaP(i)=mean(data_clean.(varname)(data_clean.(varname)>=bins(i) & data_clean.(varname)<bins(i+1)));
%     thetaN(i)=length(data_clean.(varname)(data_clean.(varname)>=bins(i) & data_clean.(varname)<bins(i+1)));
%     thetaSEM(i)=sem(data_clean.(varname)(data_clean.(varname)>=bins(i) & data_clean.(varname)<bins(i+1)));
%     insightP(i)=mean(data_clean.InsightPost(data_clean.(varname)>=bins(i) & data_clean.(varname)<bins(i+1)));
%     insightN(i)=length(data_clean.InsightPost(data_clean.(varname)>=bins(i) & data_clean.(varname)<bins(i+1)));
%     insightSEM(i)=sem(data_clean.InsightPost(data_clean.(varname)>=bins(i) & data_clean.(varname)<bins(i+1)));
% end

simpleCorPlot(data_clean.(varname),data_clean.AhaMoment)
% errorbar(thetaP,insightP,[],[],thetaP-thetaSEM,thetaP+thetaSEM)