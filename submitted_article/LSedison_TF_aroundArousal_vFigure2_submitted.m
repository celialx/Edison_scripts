%%
clear all
close all

path_fieldtrip='/Users/thandrillon/Work/local/fieldtrip/';
path_localsleep='/Users/thandrillon/WorkGit/projects/inprogress/wanderIM/localsleep';
addpath(path_fieldtrip);
addpath(path_localsleep);
ft_defaults;

path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

data_path='/Users/thandrillon/Data/LS_Edison/EDF_fixed';
files=dir([data_path filesep '*.edf']);

ColorsGroup=[93 175 117;
    133 189 181;
    67 115 128]/256;

ColorsSolver=[196 146 46 ; 177 68 22]/256;
data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).

%% Gather timing all arousals
totPerm=100;
%%
all_Files=[];
doERP=1;
all_TF_arousals=[];
StageDrop=[];
all_ERP=[];
all_maxAbs=[];
InsigthDrop=[];
nc=0;
totPerm=100;
all_TF_arousals_perm=[];
StageDrop_perm=[];
PptionSleep_after=[];
PptionSleep_after_perm=[];
PptionSleep_before=[];
PptionSleep_before_perm=[];
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
        
        ArousalFlag=T.Stade(1:end-1)~=0 & T.Stade(2:end)==0 & ~isnan(T.Stade(1:end-1));
arousals=find(ArousalFlag);
ArousalFlag=T.Stade(1:end-1)~=0 & T.Stade(2:end)==0 & ~isnan(T.Stade(1:end-1));
start_Arousal=T.TimeID(ArousalFlag)*hdr.Fs;
Beg=T.TimeID(find(~isnan(T.Stade)))*hdr.Fs; Beg=Beg(1);
diffArousals=diff([Beg ; start_Arousal])/hdr.Fs;
start_Arousal(diffArousals<60)=[];
arousals(diffArousals<60)=[];     
if isempty(start_Arousal)
    continue;
end
        %%% consider only the pause
        %     hypnogram = import_hypnogram([Folder_Name filesep File_Name(1:4) '_hypnogram.TXT']);
        addpath('../TF')
        cfg=[];
        cfg.trialfun             = 'LSedison_trialfun_arousals';
        cfg.dataset             = [Folder_Name filesep  File_Name];
        cfg.trialdef.prestim    = 2*60;
        cfg.trialdef.poststim   = .5*60;
        cfg.demean          = 'yes';
        cfg.start_Arousal   = start_Arousal;
        cfg = ft_definetrial(cfg);
        cfg.channel={'Fp1-A2','C3-A2','O1-A2'};
        data                   = ft_preprocessing(cfg); % read raw data
        
        

        keep=ones(1,length(data.trial));
        maxAbs = [];
        for k=1:length(data.trial)
            maxAbs(k,:)     = max(abs(data.trial{k}(:,data.time{1}<0)),[],2);
        end
        all_maxAbs=[all_maxAbs ; maxAbs];
        threshArt           = 750;
        nc=nc+1;
        pption_rejTr(nc,:)  = [sum(maxAbs(:,3)>threshArt) size(maxAbs,1) 100*mean(maxAbs(:,3)>threshArt)];
        fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(maxAbs(:,3)>threshArt),100*mean(maxAbs(:,3)>threshArt),threshArt)
        
        arousals(maxAbs(:,3)>threshArt)=[];
        keep(maxAbs(:,3)>threshArt)=0;
        %          trl = LSedison_trialfun_arousals(cfg)
        %         start_drop=(T.TimeID(find(T.Ball(:,1)))-1)*hdr.Fs;
        if isempty(find(keep))
            continue;
        end
        
        cfg              = [];
        cfg.trials       = find(keep);
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
        cfg.toi          = [-60:0.2:10];                         % time
        cfg.keeptrials  = 'yes';
        TFRhann = ft_freqanalysis(cfg, data);
        
        temp=10*(log10(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(:,:,:,:),3),[1 1 size(TFRhann.powspctrm,3) 1])));
        
        all_TF_arousals=cat(1,all_TF_arousals,mean(temp,1));
        all_Files=[all_Files ; repmat({SubdID},1,1)];
        if ~isempty(arousals)
            drop_indexes=T.Epoch(arousals)-1;
            temp_StageDrop=[];
            temp_PptionSleep_after=[];
            temp_PptionSleep_before=[];
            for k=1:length(arousals)
                this_stage_drop=max(unique(T.Stade(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
                this_stage_drop2=max(unique(T.N1(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
                temp_StageDrop=[temp_StageDrop ; [this_stage_drop~=0 this_stage_drop2]];
                temp_PptionSleep_after=[temp_PptionSleep_after ; [mean(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==0) mean(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==1) mean(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==2) ...
                    sum(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==0) sum(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==1) sum(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==2)]];
                
                temp_PptionSleep_before=[temp_PptionSleep_before ; [mean(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==0) mean(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==1) mean(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==2) ...
                    sum(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==0) sum(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==1) sum(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==2)]];
            end
            StageDrop=[StageDrop ; mean(temp_StageDrop,1)];
            this_line=match_str(data_clean.Sujet,SubdID);
            if isempty(this_line)
                InsigthDrop=[InsigthDrop  NaN];
            else
                InsigthDrop=[InsigthDrop  data_clean.InsightPost(this_line)];
            end
            PptionSleep_after=[PptionSleep_after ; mean(temp_PptionSleep_after,1)];
            PptionSleep_before=[PptionSleep_before ; mean(temp_PptionSleep_before(1,:),1)];
            
            Beg_Task=T.TimeID(find(T.Start(:,1)));
            End_Task=T.TimeID(find(T.End(:,1)));
            ASSM_Score_Time=T.TimeID(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
            tmprnd=round(rand(1,totPerm)*length(ASSM_Score_Time));
            while min(tmprnd)==0
                tmprnd=round(rand(1,totPerm)*length(ASSM_Score_Time));
            end
            all_rand_time=ASSM_Score_Time(tmprnd);
            %             while nPerm<totPerm
            %                 remaining_arousals=all_time_arousals(all_time_arousals(:,1)~=nF,:);
            %                 randorder=randperm(size(remaining_arousals,1));
            %                 this_rand_time=remaining_arousals(randorder(1,1),3)+Beg_Task;
            %                 if isempty(find(this_rand_time==ASSM_Score_Time))
            %                     continue;
            %                 else
            %                     nPerm=nPerm+1;
            %                 end
            %                 if (this_rand_time-60)<Beg_Task
            %                     continue;
            %                 end
            %                 all_rand_time=[all_rand_time this_rand_time];
            %             end
            cfg=[];
            cfg.trialfun             = 'LSedison_trialfun_arousals_perm';
            cfg.dataset             = [Folder_Name filesep  File_Name];
            cfg.trialdef.prestim    = 2*60;
            cfg.trialdef.poststim   = .5*60;
            cfg.demean              = 'yes';
            cfg.AllDrop             = all_rand_time;
            cfg = ft_definetrial(cfg);
            cfg.channel={'Fp1-A2','C3-A2','O1-A2'};
            data_perm               = ft_preprocessing(cfg); % read raw data
            
            keep_perm=ones(1,length(data_perm.trial));
            maxAbs = [];
            for k=1:length(data_perm.trial)
                maxAbs(k,:)     = max(abs(data_perm.trial{k}(:,data_perm.time{1}<0)),[],2);
            end
            all_rand_time(maxAbs(:,3)>threshArt)=[];
            keep_perm(maxAbs(:,3)>threshArt)=0;
            
            cfg              = [];
            cfg.trials       = find(keep_perm);
            cfg.output       = 'pow';
            cfg.channel      = 'all';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.foi          = 0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
            cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
            cfg.toi          = [-60:0.2:10];                         % time
            cfg.keeptrials  = 'yes';
            TFRhann_perm = ft_freqanalysis(cfg, data_perm);
            
            temp_perm=10*(log10(TFRhann_perm.powspctrm./repmat(mean(TFRhann_perm.powspctrm(:,:,:,:),3),[1 1 size(TFRhann_perm.powspctrm,3) 1])));
            all_TF_arousals_perm=cat(1,all_TF_arousals_perm,(nanmean(temp_perm,1)));
            
            temp_PptionSleep_after_perm=[];
            temp_PptionSleep_before_perm=[];
            temp_StageDrop_perm=[];
            for j=1:length(all_rand_time)
                this_stage_drop_perm=max(unique(T.Stade(T.Epoch==T.Epoch(all_rand_time(j))-1 | T.Epoch==T.Epoch(all_rand_time(j)))));
                this_stage_drop2_perm=max(unique(T.N1(T.Epoch==T.Epoch(all_rand_time(j))-1 | T.Epoch==T.Epoch(all_rand_time(j)))));
                temp_StageDrop_perm=[temp_StageDrop_perm ; [this_stage_drop_perm~=0 this_stage_drop2_perm]];
                
                drop_indexes=T.Epoch(all_rand_time(j))-1;
                temp_PptionSleep_after_perm=[temp_PptionSleep_after_perm ; [mean(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==0) mean(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==1) mean(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==2) ...
                    sum(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==0) sum(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==1) sum(T.Stade(T.Epoch>drop_indexes(1) & ~isnan(T.Stade))==2)]];
                
                temp_PptionSleep_before_perm=[temp_PptionSleep_before_perm ; [mean(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==0) mean(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==1) mean(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==2) ...
                    sum(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==0) sum(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==1) sum(T.Stade(T.Epoch<drop_indexes(1) & ~isnan(T.Stade))==2)]];
            end
            StageDrop_perm=[StageDrop_perm ; nanmean(temp_StageDrop_perm,1)];
            PptionSleep_after_perm=[PptionSleep_after_perm ; nanmean(temp_PptionSleep_after_perm,1)];
            PptionSleep_before_perm=[PptionSleep_before_perm ; nanmean(temp_PptionSleep_before_perm,1)];
        end
        
    end
end
%     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')

%%
Sumarousals=[];
for nF=1:length(all_Files)
    thisF=match_str(data_clean.Sujet,all_Files{nF});
    if ~isempty(thisF)
        if data_clean.InsightPre(thisF)
        SumDrop(nF)=NaN;
        else
        SumDrop(nF)=sum(~isnan(table2array(data_clean(thisF,51:54))));
        end
    else
        Sumarousals(nF)=NaN;
    end
end

%%
freqs=TFRhann.freq;
times=TFRhann.time;
% temp_toplot=squeeze(mean(all_TF_arousals(:,3,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,1);
% % caxis([-4 3])
% colorbar;
%
% xlabel('Time from Drop (s)')
% ylabel('Frequency (Hz)')
% title(sprintf('N=%g arousals',size(all_TF_arousals,1)))
% format_fig;
%%
% figure;
% subplot(1,2,1);
% temp_toplot=squeeze(mean(all_TF_arousals(StageDrop==0,3,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,0);
% % caxis([-4 4])
% colorbar;
%
% subplot(1,2,2);
% temp_toplot=squeeze(mean(all_TF_arousals(StageDrop~=0,3,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,0);
% % caxis([-4 4])
% colorbar;

%%
% figure
% temp_toplot=squeeze(mean(all_TF_arousals(~isnan(Sumarousals),COI,freqs>=FOI(1) & freqs<=FOI(2),times>-50),3))-...
%     squeeze(mean(all_TF_arousals_perm(~isnan(Sumarousals),COI,freqs>=FOI(1) & freqs<=FOI(2),times>-50),3));
% simpleTplot(times(times>-50),temp_toplot,0,'k',[2 0.05 0.05 1000],'-',0.5,1,10,1,1);


%%
Freqs=[1 4; 4 8; 8 12];
ylims=[-1.4 1.2; -4 0; -2 1];
TitleFreqs={'\delta','\theta','\alpha','\sigma'};

COI=3;
figure; set(gcf,'Position',[57   338   400   349]);
% for nF=1:size(Freqs,1)
%     subplot(1,3,nF);
nF=1;
FOI=Freqs(nF,:);
temp_toplot=squeeze(mean(all_TF_arousals(~isnan(Sumarousals),COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,10,1,1);

temp_toplot=squeeze(mean(all_TF_arousals_perm(~isnan(Sumarousals),COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
plot(times,mean(temp_toplot,1),'Color',[1 1 1]*0.4,'LineStyle','-','LineWidth',3);
plot(times,mean(temp_toplot,1)+sem(temp_toplot,1),'Color',[1 1 1]*0.8,'LineStyle','--','LineWidth',2);
plot(times,mean(temp_toplot,1)-sem(temp_toplot,1),'Color',[1 1 1]*0.8,'LineStyle','--','LineWidth',2);

% temp_toplot=squeeze(mean(all_TF_arousals(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
% simpleTplot(times,mean(temp_toplot),0,ColorsGroup(1,:),0,'-',0.5,1,10,1,3);
%
% temp_toplot=squeeze(mean(all_TF_arousals(StageDrop==1,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
% simpleTplot(times,mean(temp_toplot),0,ColorsGroup(3,:),0,'-',0.5,1,10,1,3);
% simpleTplot(times,temp_toplot,0,ColorsGroup(3,:),0,'-',0.5,1,10,1,1);

xlim([-50 5])
% ylim(ylims(nF,:))
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
% title(TitleFreqs{nF})
line([0 0],ylims,'Color',[1 1 1]*.7,'LineStyle','--');
set(gca,'LineWidth',2);
% end
%  hold on;
%  temp_toplot=squeeze(mean(all_TF_arousals(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,0,1,1);
%  temp_toplot=squeeze(mean(all_TF_arousals(StageDrop~=0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'b',0,'-',0.5,1,0,1,1);
%
% line([-19.801 -0.602],[1 1]*ylims(1,2),'LineWidth',2,'Color','k');
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_TimePlot_v2.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_TimePlot_v2.png'],'-r 300')



%% Linear fit
COI=3;
begWinSlop=-50;
FOI=[1 4]; %Freqs(1,:);
vec_power=[];
vec_time=[];
coeff_polyfit{1}=[];
coeff_polyfit{2}=[];
temp_toplot=squeeze(mean(all_TF_arousals(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
temp_toplot_perm=squeeze(mean(all_TF_arousals_perm(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
for k=1:size(temp_toplot,1)
    vec_power=[vec_power temp_toplot(k,times>begWinSlop & times<-2)];
    vec_time=[vec_time times(times>begWinSlop & times<-2)];
    [p,S,mu] = polyfit(times(times>begWinSlop & times<-2),temp_toplot(k,times>begWinSlop & times<-2),1);
    
    coeff_polyfit{1}=[coeff_polyfit{1} ; p];
    
    
    [p,S,mu] = polyfit(times(times>begWinSlop & times<-2),temp_toplot_perm(k,times>begWinSlop & times<-2),1);
    coeff_polyfit{2}=[coeff_polyfit{2} ; p];
end

figure;
set(gcf,'Position',[ 287   325   299   419]);
pV_slope=[];
format_fig;
StageDrop2=StageDrop(:,2);
for nCond=1:2
    temp_toplot=coeff_polyfit{nCond}(~isnan(Sumarousals),1);
    [pV_slope(nCond)]=signrank(temp_toplot,0);
    %%%%% Plot
    hold on;
    Pos=nCond; data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
    %         colorBar=[ColorsGroup(2*nCond-1,:); 0 0 0];
    if nCond==1
        colorBar=[0.2 .2 0.2; 0 0 0];
    else
        colorBar=[1 1 1; 0 0 0];
    end
    line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
    line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
    
    patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
        [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
    
    if nCond==1
        xspread=(rand(1,length(data(StageDrop2(~isnan(Sumarousals))==0)))-0.5)*widthBar/6+nCond+2.5*0.1*widthBar;
        scatter(xspread,data(StageDrop2(~isnan(Sumarousals))==0),'Marker',markerType,'MarkerFaceColor',ColorsGroup(1,:),'MarkerEdgeColor',ColorsGroup(1,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
        
        
        xspread=(rand(1,length(data(StageDrop2(~isnan(Sumarousals))~=0)))-0.5)*widthBar/6+nCond+0.2+2.5*0.1*widthBar;
        scatter(xspread,data(StageDrop2(~isnan(Sumarousals))~=0),'Marker',markerType,'MarkerFaceColor',ColorsGroup(3,:),'MarkerEdgeColor',ColorsGroup(3,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    else
        xspread=(rand(1,length(data))-0.5)*widthBar/6+nCond+2.5*0.1*widthBar;
        scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
end
xlim([.5 2.8])
ylim([-3     2.2])
set(gca,'XTick',1:2,'XTickLabel',{'Real','Perm'},'LineWidth',2)
ylabel('Linear Fit - Slope Delta Before Drop')
line(xlim,[0 0],'Color','k','LineStyle','--');

[pV_slope(3)]=signrank(coeff_polyfit{1}((~isnan(Sumarousals)),1),coeff_polyfit{2}((~isnan(Sumarousals)),1));

% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_Slope_v2.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_Slope_v2.png'],'-r 300')

% A=[coeff_polyfit{1}(:,1) ; coeff_polyfit{2}(:,1)];
% B=[[zeros(size(coeff_polyfit{1},1),1) ; ones(size(coeff_polyfit{1},1),1)] [StageDrop2 ; nan(size(coeff_polyfit{1},1),1)]];
% anova1(A,B(:,1));
%%
figure;
set(gcf,'Position',[ 287   325   299   419]);
pV_duration=[];
format_fig;
StageDrop2=StageDrop(:,2);
for nCond=1:2
    if nCond==1
        temp_toplot=100*sum(PptionSleep_after((~isnan(Sumarousals)),2:3),2);
    else
        temp_toplot=100*sum(PptionSleep_after_perm((~isnan(Sumarousals)),2:3),2);
    end
    [pV_post(nCond)]=signrank(temp_toplot,0);
    %%%%% Plot
    hold on;
    Pos=nCond; data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
    %         colorBar=[ColorsGroup(2*nCond-1,:); 0 0 0];
    if nCond==1
        colorBar=[0.4 .4 0.4; 0 0 0];
    else
        colorBar=[1 1 1; 0 0 0];
    end
    line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
    line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
    line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
    
    patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
        [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
    
    if nCond==1
        xspread=(rand(1,length(data(StageDrop2(~isnan(Sumarousals))==0)))-0.5)*widthBar/6+nCond+2.5*0.1*widthBar;
        scatter(xspread,data(StageDrop2(~isnan(Sumarousals))==0),'Marker',markerType,'MarkerFaceColor',ColorsGroup(1,:),'MarkerEdgeColor',ColorsGroup(1,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
        
        
        xspread=(rand(1,length(data(StageDrop2(~isnan(Sumarousals))~=0)))-0.5)*widthBar/6+nCond+0.2+2.5*0.1*widthBar;
        scatter(xspread,data(StageDrop2(~isnan(Sumarousals))~=0),'Marker',markerType,'MarkerFaceColor',ColorsGroup(3,:),'MarkerEdgeColor',ColorsGroup(3,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    else
        xspread=(rand(1,length(data))-0.5)*widthBar/6+nCond+2.5*0.1*widthBar;
        scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
end
    [pV_post(3)]=signrank(sum(PptionSleep_after((~isnan(Sumarousals)),2:3),2),sum(PptionSleep_after_perm((~isnan(Sumarousals)),2:3),2));

    xlim([.5 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Real','Perm'},'LineWidth',2)
ylabel('Ption Sleep Post Drop')
line(xlim,[0 0],'Color','k','LineStyle','--');

% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_PostDrop_SleepPption.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_PostDrop_SleepPption.png'],'-r 300')

A2=[sum(PptionSleep_after(:,2:3),2) ; sum(PptionSleep_after_perm(:,2:3),2)];
B2=[[zeros(size(coeff_polyfit{1},1),1) ; ones(size(coeff_polyfit{1},1),1)] [StageDrop2 ; nan(size(coeff_polyfit{1},1),1)]];

%%
figure;
set(gcf,'Position',[ 287   325   299   419]);
for nCond=1:2
    temp_toplot=sum(PptionSleep_after(:,nCond+4),2);
    %%%%% Plot
    for k=1:2
        hold on;
        Pos=nCond+(k-2)*0.2; data=temp_toplot(InsigthDrop==k-1); widthLine=2; widthBar=0.8; sizeDot=400; markerType='o';
        %         colorBar=[ColorsGroup(2*nCond-1,:); 0 0 0];
        if nCond==1
            colorBar=[ColorsSolver(k,:); 0 0 0];
        else
            colorBar=[ColorsSolver(k,:); 0 0 0];
        end
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
        line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        
        patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
            [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
        
        xspread=(rand(1,length(data))-0.5)*widthBar/6+nCond+2.5*0.12*widthBar;
        scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
end
xlim([.5 2.5])
set(gca,'XTick',1:2,'XTickLabel',{'N1','N2'},'LineWidth',2)
ylabel('Ption Post Drop')
line(xlim,[0 0],'Color','k','LineStyle','--');
format_fig;

% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_PostDrop_SleepPption_byInsight.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_PostDrop_SleepPption_byInsight.png'],'-r 300')

[pV_Ins(1)]=ranksum(sum(PptionSleep_after(InsigthDrop==0,2),2),sum(PptionSleep_after(InsigthDrop==1,2),2));
[pV_Ins(2)]=ranksum(sum(PptionSleep_after(InsigthDrop==0,2),3),sum(PptionSleep_after(InsigthDrop==1,3),2));


%%
figure;
set(gcf,'Position',[ 287   325   299   419]);
for nCond=1:2
    temp_toplot=sum(PptionSleep_before(:,nCond+4),2)+sum(PptionSleep_after(:,nCond+4),2);
    %%%%% Plot
    for k=1:2
        all_toplots{nCond,k}=temp_toplot(InsigthDrop==k-1);
        hold on;
        Pos=nCond+(k-2)*0.2; data=temp_toplot(InsigthDrop==k-1); widthLine=2; widthBar=0.8; sizeDot=400; markerType='o';
        %         colorBar=[ColorsGroup(2*nCond-1,:); 0 0 0];
        if nCond==1
            colorBar=[ColorsSolver(k,:); 0 0 0];
        else
            colorBar=[ColorsSolver(k,:); 0 0 0];
        end
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
        line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        
        patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
            [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
        
        xspread=(rand(1,length(data))-0.5)*widthBar/6+nCond+2.5*0.12*widthBar;
        scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
end
xlim([.5 2.5])
set(gca,'XTick',1:2,'XTickLabel',{'N1','N2'})
ylabel('Ption Post Drop')
line(xlim,[0 0],'Color','k','LineStyle','--');
format_fig;


[pV_Ins2(1)]=ranksum(all_toplots{1,1},all_toplots{1,2});
[pV_Ins2(2)]=ranksum(all_toplots{2,1},all_toplots{2,2});


%%
figure;
set(gcf,'Position',[ 287   325   350   419]);
hold on;
bar(1,36.67,'BarWidth',0.65,'FaceColor','k','EdgeColor','k','LineWidth',3);
bar(2,61.7,'BarWidth',0.65,'FaceColor',[1 1 1]*0.7,'EdgeColor','k','LineWidth',3);
bar(3,22.5,'BarWidth',0.65,'FaceColor',[1 1 1]*1,'EdgeColor','k','LineWidth',3);

xlim([.2 3.8])
ylim([0 80])
set(gca,'XTick',1:3,'XTickLabel',{'N1','Drop','ø Drop'},'XColor','k','YColor','k','LineWidth',2)
ylabel('Hypnagogia (%)')
line([1.5 1.5],ylim,'Color','k','LineStyle',':','LineWidth',2);
format_fig;

% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_Hypngogia.eps'],'-r 300')
% export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_Hypngogia.png'],'-r 300')

%%
%%
Freqs=[1 4; 4 8; 8 12];
ylims=[-1.4 1.2; -4 0; -2 1];
TitleFreqs={'\delta','\theta','\alpha','\sigma'};

COI=3;
figure; set(gcf,'Position',[57   338   400   349]);
% for nF=1:size(Freqs,1)
%     subplot(1,3,nF);
nF=1;
FOI=Freqs(nF,:);
temp_toplot=squeeze(mean(all_TF_arousals(InsigthDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,10,1,1);

temp_toplot=squeeze(mean(all_TF_arousals(InsigthDrop==1,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,10,1,1);

xlim([-50 5])
% ylim(ylims(nF,:))
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
% title(TitleFreqs{nF})
line([0 0],ylims,'Color',[1 1 1]*.7,'LineStyle','--');
set(gca,'LineWidth',2);
