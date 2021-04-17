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
files=dir([data_path filesep '*.edf']);

ColorsGroup=[55 179 111;
    115 191 181;
    50 116 130]/256;

ColorsSolver=[206 144 45 ; 192 50 2]/256;
data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).


%%
doERP=1;
all_TF_drops=[];
StageDrop=[];
all_ERP=[];
all_maxAbs=[];
InsigthDrop=[];
nc=0;
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
        drops=find(T.Ball(:,1));
        
        keep=ones(1,length(data.trial));
        maxAbs = [];
        for k=1:length(data.trial)
            maxAbs(k,:)     = max(abs(data.trial{k}(:,data.time{1}<0)),[],2);
        end
        all_maxAbs=[all_maxAbs ; maxAbs];
        threshArt           = 150;
        nc=nc+1;
        pption_rejTr(nc,:)  = [sum(maxAbs(:,3)>threshArt) size(maxAbs,1) 100*mean(maxAbs(:,3)>threshArt)];
        fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(maxAbs(:,3)>threshArt),100*mean(maxAbs(:,3)>threshArt),threshArt)
        
        drops(maxAbs(:,3)>threshArt)=[];
        %          trl = LSedison_trialfun_drops(cfg)
        %         start_drop=(T.TimeID(find(T.Ball(:,1)))-1)*hdr.Fs;
        
        
        cfg              = [];
        cfg.trials       = find(keep);
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 0.5:0.2:30;                         % analysis 2 to 30 Hz in steps of .2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*6;   % length of time window = 0.5 sec
        cfg.toi          = [-50:0.2:10];                         % time
        cfg.keeptrials  = 'yes';
        TFRhann = ft_freqanalysis(cfg, data);
        
        temp=10*(log10(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(:,:,:,:),3),[1 1 size(TFRhann.powspctrm,3) 1])));
        
        all_TF_drops=cat(1,all_TF_drops,temp);
        
        if ~isempty(drops)
            drop_indexes=T.Epoch(drops)-1;
            for k=1:length(drops)
                this_stage_drop=max(unique(T.Stade(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
                StageDrop=[StageDrop this_stage_drop];
                this_line=match_str(data_clean.Sujet,SubdID);
                if isempty(this_line)
                    InsigthDrop=[InsigthDrop  NaN];
                else
                    InsigthDrop=[InsigthDrop  data_clean.InsightPost(this_line)];
                end
            end
        end
    end
    %     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
end

%%
freqs=TFRhann.freq;
times=TFRhann.time;
% temp_toplot=squeeze(mean(all_TF_drops(:,3,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,1);
% % caxis([-4 3])
% colorbar;
% 
% xlabel('Time from Drop (s)')
% ylabel('Frequency (Hz)')
% title(sprintf('N=%g drops',size(all_TF_drops,1)))
% format_fig;
%%
% figure;
% subplot(1,2,1);
% temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,3,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,0);
% % caxis([-4 4])
% colorbar;
% 
% subplot(1,2,2);
% temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,3,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,0);
% % caxis([-4 4])
% colorbar;

%%
Freqs=[1 4; 4 8; 8 12];
ylims=[-0.55 1.65; -4 0; -2 1];
TitleFreqs={'\delta','\theta','\alpha','\sigma'};

COI=2;
figure; set(gcf,'Position',[57   338   650   349]);
% for nF=1:size(Freqs,1)
%     subplot(1,3,nF);
nF=1;
FOI=Freqs(nF,:);
temp_toplot=squeeze(mean(all_TF_drops(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,10,1,1);

% temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
% simpleTplot(times,mean(temp_toplot),0,ColorsGroup(1,:),0,'-',0.5,1,10,1,3);
% 
% temp_toplot=squeeze(mean(all_TF_drops(StageDrop==1,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
% simpleTplot(times,mean(temp_toplot),0,ColorsGroup(3,:),0,'-',0.5,1,10,1,3);
% simpleTplot(times,temp_toplot,0,ColorsGroup(3,:),0,'-',0.5,1,10,1,1);

xlim([-50 5])
ylim(ylims(nF,:))
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
% title(TitleFreqs{nF})
line([0 0],ylims,'Color',[1 1 1]*.7,'LineStyle','--');

% end
%  hold on;
%  temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,0,1,1);
%  temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'b',0,'-',0.5,1,0,1,1);
%

export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_TimePlot.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_TimePlot.png'],'-r 300')

%% Linear fit
COI=2;
for nF=1 %:size(Freqs,1)
    FOI=Freqs(nF,:);
    vec_power=[];
    vec_time=[];
    coeff_polyfit{nF}=[];
    temp_toplot=squeeze(mean(all_TF_drops(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
    for k=1:size(temp_toplot,1)
        vec_power=[vec_power temp_toplot(k,times>-116 & times<-2)];
        vec_time=[vec_time times(times>-116 & times<-2)];
        [p,S,mu] = polyfit(times(times>-116 & times<-2),temp_toplot(k,times>-116 & times<-2),1);
        
        coeff_polyfit{nF}=[coeff_polyfit{nF} ; p];
    end
end

figure;
set(gcf,'Position',[ 287   325   299   419]);
pV_slope=[];
format_fig;
for nF=1 %:size(Freqs,1)
    for nCond=1:2
        temp_toplot=coeff_polyfit{nF}(StageDrop==nCond-1,1);
        [pV_slope(nF)]=signrank(temp_toplot,0);
        %%%%% Plot
        hold on;
        Pos=nCond; data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
        colorBar=[ColorsGroup(2*nCond-1,:); 0 0 0];
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
        line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        
        patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
            [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
        
        xspread=(rand(1,length(data))-0.5)*widthBar/6+nCond+2.5*0.1*widthBar;
        scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
    end
end
xlim([.5 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Wake','Sleep'})
ylabel('Linear Fit - Slope Delta Before Drop')
line(xlim,[0 0],'Color','k','LineStyle','--');


export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_Slope.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_AroundDrop_DeltaPower_Slope.png'],'-r 300')

%%
cmap=[242,114,8;22,92,125]/255;
FOI=[1 4];
COI=3;
figure;
temp_toplot=squeeze(mean(all_TF_drops(InsigthDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,ColorsSolver(2,:),0,'-',0.5,1,10,1,1);

temp_toplot=squeeze(mean(all_TF_drops(InsigthDrop==1,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,ColorsSolver(1,:),0,'-',0.5,1,10,1,1);

xlim([-50 10])
ylim([-1.6 2.8])
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
title(TitleFreqs{nF})
line([0 0],ylims,'Color',[1 1 1]*.7,'LineStyle','--');

