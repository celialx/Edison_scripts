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
        for k=1:length(data.trial)
            all_ERP=[all_ERP ; data.trial{k}(3,:)];
            
            if max(abs(data.trial{k}(3,data.time{k}<0)))>100
                keep(k)=0;
            end
        end
        if isempty(find(keep))
            continue;
        end
        drops(keep==0)=[];
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
        cfg.toi          = [-2*60:2:0.5*60];                         % time
        cfg.keeptrials  = 'yes';
        TFRhann = ft_freqanalysis(cfg, data);
        
        temp=(log(TFRhann.powspctrm./repmat(mean(TFRhann.powspctrm(:,:,:,:),3),[1 1 size(TFRhann.powspctrm,3) 1])));
        
        all_TF_drops=cat(1,all_TF_drops,temp);
        
        if ~isempty(drops)
            drop_indexes=T.Epoch(drops)-1;
            for k=1:length(drops)
                this_stage_drop=max(unique(T.Stade(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
                StageDrop=[StageDrop this_stage_drop];
            end
        end
    end
    %     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
end

%%
freqs=TFRhann.freq;
times=TFRhann.time;
temp_toplot=squeeze(mean(all_TF_drops(:,3,:,:),1));
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
temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,3,:,:),1));
h=simpleTFplot(temp_toplot,freqs,times,0,0);
caxis([-4 4])

subplot(1,2,2);
temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,3,:,:),1));
h=simpleTFplot(temp_toplot,freqs,times,0,0);
caxis([-4 4])

%%
Freqs=[1 4; 4 8; 8 12 ; 15 20];
TitleFreqs={'detla','theta','alpha','sigma'};

COI=2;
figure;
for nF=1:size(Freqs,1)
    subplot(1,4,nF);
    FOI=Freqs(nF,:);
    temp_toplot=squeeze(mean(all_TF_drops(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%     temp_toplot=temp_toplot-repmat(temp_toplot(:,times==-5),1,size(temp_toplot,2)); %nanzscore(temp_toplot,[],2);
    simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,0,1,1);

xlim([-120 10])
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
title(sprintf('N=%g drops',size(all_TF_drops,1)))
end
%  hold on;
%  temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,0,1,1);
%  temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
%  simpleTplot(times,temp_toplot,0,'b',0,'-',0.5,1,0,1,1);
%

%% Linear fit
COI=3;
for nF=1:size(Freqs,1)
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
set(gcf,'Position',[287   325   573   419]);
format_fig;
for nF=1:size(Freqs,1)
        temp_toplot=coeff_polyfit{nF}(:,1);
        [pV_slope(nF)]=signrank(temp_toplot,0);
        %%%%% Plot
        hold on;
        Pos=nF; data=temp_toplot; widthLine=2; widthBar=1.2; sizeDot=400; markerType='o';
        colorBar=[[1 1 1]*0.7 ; [1 1 1]*0.7];
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
        line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
        
        patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
            [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');
        
        xspread=(rand(1,length(data))-0.5)*widthBar/6+nF+2.5*0.1*widthBar;
        scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);

end
xlim([.5 4.8])
set(gca,'XTick',1:4,'XTickLabel',{'\delta','\theta','\alpha','\sigma'})
ylabel('Linear Fit - Slope')
line(xlim,[0 0],'Color','k','LineStyle','--');
%%
cmap=[242,114,8;22,92,125]/255;
FOI=[1 4];
COI=2;
figure;
temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,cmap(1,:),0,'-',0.5,1,0,1,1);

temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,cmap(2,:),0,'-',0.5,1,0,1,1);

xlim([-60 20])
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
title(sprintf('N=%g drops',size(all_TF_drops,1)))