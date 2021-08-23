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
all_fooof_drops=cell(3,2);
all_spec_drops=cell(3,1);
StageDrop=[];
all_hypno_ASSM=[];
all_hypno_Stages=[];
all_hypnoD_Stages=[];
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
    if exist([Folder_Name filesep '..' filesep 'StanfordStages' filesep File_Name(1:end-4) '.hypnogram.txt'])==0
        warning('Stanford-stages file missing');
        continue;
    end
    fprintf('... processing %s (%g/%g)\n',File_Name,nF,length(files));
    load([Folder_Name filesep 'T_' SubdID '.mat'])
    hypnodensity = import_stage_hypnodensity([Folder_Name filesep '..' filesep 'StanfordStages' filesep File_Name(1:end-4) '.hypnodensity.txt']);
    hypnogram = import_stage_hypnogram([Folder_Name filesep '..' filesep 'StanfordStages' filesep File_Name(1:end-4) '.hypnogram.txt']);
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([Folder_Name filesep File_Name]);
    
    %     figure;
    ASSM_Score=T.Stade;
    Beg_Task=T.TimeID(find(T.Start(:,1)));
    End_Task=T.TimeID(find(T.End(:,1)));
    ASSM_Score=ASSM_Score(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    ASSM_Score_Time=T.TimeID(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    
    
    hypnogram_perSec=reshape(repmat(hypnogram,1,15)',1,numel(repmat(hypnogram,1,15)));
    time_hypno_perSec=1:length(hypnogram_perSec);
    Stanford_Score=hypnogram_perSec(time_hypno_perSec>=Beg_Task & time_hypno_perSec<=End_Task);
    Stanford_Score_Time=time_hypno_perSec(time_hypno_perSec>=Beg_Task & time_hypno_perSec<=End_Task);
    Stanford_HypnoD=nan(5,size(Stanford_Score,2));
    for j=1:size(hypnodensity,2)
        temparray=table2array(hypnodensity(:,j));
        temparray2=reshape(repmat(temparray,1,15)',1,numel(repmat(temparray,1,15)));
        Stanford_HypnoD(j,:)=temparray2(time_hypno_perSec>=Beg_Task & time_hypno_perSec<=End_Task);
    end
    %     scores=find(T.StadeL~='?');
    drops=find(T.Ball(:,1));
    if ~isempty(drops)
        drop_indexes=T.Epoch(find(T.Ball(:,1)))-1;
        
        for k=1:length(drops)
            time_drop=T.TimeID(drops(k));
            
            this_stage_drop=max(unique(T.Stade(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
            
            StageDrop=[StageDrop this_stage_drop];
            temp_hypno_ASSM=ASSM_Score(find(ASSM_Score_Time==time_drop)+(-120:30));
            temp_hypno_Stages=Stanford_Score(find(Stanford_Score_Time==time_drop)+(-120:30));
            temp_hypnoD_Stages=Stanford_HypnoD(:,find(Stanford_Score_Time==time_drop)+(-120:30));
            
            all_hypno_ASSM=[all_hypno_ASSM ; temp_hypno_ASSM'];
            all_hypno_Stages=[all_hypno_Stages ; temp_hypno_Stages];
            all_hypnoD_Stages=cat(3,all_hypnoD_Stages, temp_hypnoD_Stages);
        end
    end
    nc=nc+1;
    perS_drops(nc)=length(drops);
    perS_deepestSt(nc)=max(ASSM_Score);
    
end

%%
cmap=[242,114,8;22,92,125]/255;
figure;
simpleBarPlot(0.8,mean(perS_drops~=0 & perS_deepestSt==0),cmap(1,:),0.35,'k',[],3);
simpleBarPlot(1.2,mean(perS_drops==0 & perS_deepestSt==0),[1 1 1;cmap(1,:)],0.35,'k',[],3);

simpleBarPlot(1.8,mean(perS_drops~=0 & perS_deepestSt~=0),cmap(2,:),0.35,'k',[],3);
simpleBarPlot(2.2,mean(perS_drops==0 & perS_deepestSt~=0),[1 1 1;cmap(2,:)],0.35,'k',[],3);

ylabel('Proportion of participants')
set(gca,'XTick',1:2,'XTickLabel',{'No','Yes'});
xlabel('Slept?')
format_fig;
xlim([0.2 2.8])
%% ASSM

figure; set(gcf,'Position',[449   427   791   378]);
format_fig;
times=-120:30;
temp_plot=all_hypno_ASSM(StageDrop==0,:)==0;
% plot(times,temp_plot','k');
simpleTplot(times,temp_plot,0,cmap(1,:),0,'-',0.5,1,0,1,1);
hold on
temp_plot=all_hypno_ASSM(StageDrop~=0,:)==0;
simpleTplot(times,temp_plot,0,cmap(2,:),0,'-',0.5,1,0,1,1);
xlim([-120 20])
ylim([0 1])
set(gca,'YTick',0:0.2:1)
ylabel('Probability Wake')
xlabel('s from drop')
line([0 0],ylim,'Color','k','LineStyle',':')
title('AASM')

%% Stages
cmap=[242,114,8;22,92,125]/255;
maxSize=288;
figure; set(gcf,'Position',[449   427   791   378]);
hold on;
for nx=1:151
    temp=all_hypno_ASSM(StageDrop==0,nx);
    for j=0:2
        if mean(temp==j)==0
            continue;
        end
    scatter(times(nx),j-0.05,'SizeData',maxSize*mean(temp==j),'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerFaceAlpha',0.5);
    end
    
        temp=all_hypno_ASSM(StageDrop~=0,nx);
    for j=0:2
        if mean(temp==j)==0
            continue;
        end
    scatter(times(nx),j+0.05,'SizeData',maxSize*mean(temp==j),'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerFaceAlpha',0.5);
    end
end
ylim([-0.12 2.12])

xlim([-120 20])
ylabel('Sleep Stage')
xlabel('s from drop')
title('Stages')
set(gca,'YTick',0:2,'YTickLabel',{'WK','N1','N2'});

line(xlim,[1 1]*0.5,'Color','k','LineStyle','--')
line([0 0],ylim,'Color','k','LineStyle',':')
format_fig;
%%
figure; set(gcf,'Position',[449   427   791   378]);
format_fig;
times=-120:30;
temp_plot=squeeze(all_hypnoD_Stages(1,:,StageDrop==0))';
% plot(times,temp_plot','k');
simpleTplot(times,temp_plot,0,cmap(1,:),0,'-',0.5,1,0,1,1);
hold on
temp_plot=squeeze(all_hypnoD_Stages(1,:,StageDrop~=0))';
simpleTplot(times,temp_plot,0,cmap(2,:),0,'-',0.5,1,0,1,1);
xlim([-120 20])
ylim([0 1])
set(gca,'YTick',0:0.2:1)
ylabel('Probability Wake')
xlabel('s from drop')
title('Stages')


line(xlim,[1 1]*0.5,'Color','k','LineStyle','--')
line([0 0],ylim,'Color','k','LineStyle',':')


%% Stages
figure;
format_fig;
times=-120:30;
Colors=[0.8 0 0; 0 0 0.6; 0 0 0.8; 0 0 1];
for k=1:2
    if k==2
        temp_plot=squeeze(sum(all_hypnoD_Stages(k:3,:,:),1))';
    else
        temp_plot=squeeze(all_hypnoD_Stages(k,:,:))';
    end
    [~,hp(k)]=simpleTplot(times,temp_plot,0,Colors(k,:),0,'-',0.3,1,0,1,1);
end
xlim([-120 20])
title('Stages')
legend(hp,{'W','N1+N2'});
xlabel('Time (s)');
ylabel('Proability')