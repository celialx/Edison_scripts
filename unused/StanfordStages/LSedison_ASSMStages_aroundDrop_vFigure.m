%%
clear all
close all

path_raincloud='/Users/tand0009/WorkGit/projects/ext/RainCloudPlots/';
addpath(genpath(path_raincloud));

path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

save_path='/Users/tand0009/Data/LS_Edison/LocalSleep';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
data_path='/Users/tand0009/Data/LS_Edison/EDF_fixed';
files=dir([data_path filesep '*.edf']);

data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).

%% Extract all drop timings
all_time_drops=[];
all_Transitions=[];
totPerm=100;
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
     fprintf('... processing %s (%g/%g)\n',File_Name,nF,length(files));
   load([Folder_Name filesep 'T_' SubdID '.mat'])
    %     figure;
    ASSM_Score=T.Stade;
    Beg_Task=T.TimeID(find(T.Start(:,1)));
    End_Task=T.TimeID(find(T.End(:,1)));
    ASSM_Score=ASSM_Score(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    ASSM_Score_Time=T.TimeID(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    
    this_line=match_str(data_clean.Sujet,SubdID);
    if isempty(this_line)
        Insigth=NaN;
        SleepG=NaN;
    else
        Insigth=data_clean.InsightPost(this_line);
        SleepG=data_clean.SleepEdison(this_line);
    end
    drops=find(T.Ball(:,1));
    
    all_Transitions=[all_Transitions ; [Insigth SleepG ~isempty(drops) sum(diff(ASSM_Score)~=0) sum(diff(ASSM_Score)~=0 & ASSM_Score(2:end)==1) sum(diff(ASSM_Score)~=0 & ASSM_Score(2:end)==2)]];
    if ~isempty(drops)
        drop_indexes=T.Epoch(find(T.Ball(:,1)))-1;
        
        for k=1:length(drops)
            time_drop=T.TimeID(drops(k));
            
            all_time_drops=[all_time_drops ; [nF k time_drop-Beg_Task+1]];
        end
    end
end
%%
StageDrop=[];
BeforeAfterDrop=[];
BeforeAfterDrop_perm=[];

BeforeAfterDrop2=[];
BeforeAfterDrop_perm2=[];

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
    fprintf('... processing %s (%g/%g)\n',File_Name,nF,length(files));
    load([Folder_Name filesep 'T_' SubdID '.mat'])
    
    %     figure;
    ASSM_Score=T.Stade;
    Beg_Task=T.TimeID(find(T.Start(:,1)));
    End_Task=T.TimeID(find(T.End(:,1)));
    ASSM_Score=ASSM_Score(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    ASSM_Score_Time=T.TimeID(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    
    drops=find(T.Ball(:,1));
    if ~isempty(drops)
        drop_indexes=T.Epoch(find(T.Ball(:,1)))-1;
        
        for k=1 %:length(drops)
            time_drop=T.TimeID(drops(k));
            this_stage_drop=ASSM_Score(find(ASSM_Score_Time==time_drop)); %max(unique(T.Stade(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
            
            StageDrop=[StageDrop this_stage_drop];
            
            beforeDrop=[];
            afterDrop=[];
            for j=0:2
                beforeDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time<time_drop))==j);
                afterDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time>time_drop))==j);
            end
            BeforeAfterDrop=[BeforeAfterDrop ; [nF k NaN this_stage_drop beforeDrop afterDrop]];
            
            nPerm=0;
            fprintf('%2.0f%%\n',0)
            BeforeAfterDrop_temp=[];
            while nPerm<totPerm
                remaining_drops=all_time_drops(all_time_drops(:,1)~=nF,:);
                randorder=randperm(size(remaining_drops,1));
                this_rand_time=remaining_drops(randorder(1,1),3)+Beg_Task;
                if isempty(find(this_rand_time==ASSM_Score_Time))
                    continue;
                else
                    nPerm=nPerm+1;
                end
                fprintf('\b\b\b\b%2.0f%%\n',round(100*nPerm/totPerm))
                this_stage_drop=ASSM_Score(find(ASSM_Score_Time==this_rand_time));
                
                beforeDrop=[];
                afterDrop=[];
                for j=0:2
                    beforeDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time<this_rand_time))==j);
                    afterDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time>this_rand_time))==j);
                end
                BeforeAfterDrop_temp=[BeforeAfterDrop_temp ; [nF k nPerm this_stage_drop beforeDrop afterDrop]];
            end
            BeforeAfterDrop_perm=[BeforeAfterDrop_perm ; nanmean(BeforeAfterDrop_temp,1)];
            
            
            %%%%%
            beforeDrop=[];
            afterDrop=[];
            for j=0:2
                beforeDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time<time_drop & ASSM_Score_Time>time_drop-5*60))==j);
                afterDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time>time_drop & ASSM_Score_Time<time_drop+5*60))==j);
            end
            BeforeAfterDrop2=[BeforeAfterDrop2 ; [nF k NaN this_stage_drop beforeDrop afterDrop]];
            
            nPerm=0;
            fprintf('%2.0f%%\n',0)
            BeforeAfterDrop_temp=[];
            while nPerm<totPerm
                remaining_drops=all_time_drops(all_time_drops(:,1)~=nF,:);
                randorder=randperm(size(remaining_drops,1));
                this_rand_time=remaining_drops(randorder(1,1),3)+Beg_Task;
                if isempty(find(this_rand_time==ASSM_Score_Time))
                    continue;
                else
                    nPerm=nPerm+1;
                end
                fprintf('\b\b\b\b%2.0f%%\n',round(100*nPerm/totPerm))
                this_stage_drop=ASSM_Score(find(ASSM_Score_Time==this_rand_time));
                
                beforeDrop=[];
                afterDrop=[];
                for j=0:2
                    beforeDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time<this_rand_time & ASSM_Score_Time>time_drop-5*60))==j);
                    afterDrop(j+1)=mean(ASSM_Score(find(ASSM_Score_Time>this_rand_time & ASSM_Score_Time<time_drop+5*60))==j);
                end
                BeforeAfterDrop_temp=[BeforeAfterDrop_temp ; [nF k nPerm this_stage_drop beforeDrop afterDrop]];
            end
            BeforeAfterDrop_perm2=[BeforeAfterDrop_perm2 ; nanmean(BeforeAfterDrop_temp,1)];
            
        end
    end
end


%%
forANOVA=[];
forANOVA=[forANOVA ; [BeforeAfterDrop2(:,7) ones(size(BeforeAfterDrop2,1),1) ones(size(BeforeAfterDrop2,1),1) ]];
forANOVA=[forANOVA ; [BeforeAfterDrop2(:,10) 2*ones(size(BeforeAfterDrop2,1),1) ones(size(BeforeAfterDrop2,1),1) ]];
forANOVA=[forANOVA ; [BeforeAfterDrop_perm2(:,7) ones(size(BeforeAfterDrop_perm2,1),1) zeros(size(BeforeAfterDrop_perm2,1),1) ]];
forANOVA=[forANOVA ; [BeforeAfterDrop_perm2(:,10) 2*ones(size(BeforeAfterDrop_perm2,1),1) zeros(size(BeforeAfterDrop_perm2,1),1) ]];

%%
BeforeAfterDrop(:,11)=BeforeAfterDrop(:,7)+BeforeAfterDrop(:,6);
BeforeAfterDrop(:,12)=BeforeAfterDrop(:,10)+BeforeAfterDrop(:,9);

BeforeAfterDrop_perm(:,11)=BeforeAfterDrop_perm(:,7)+BeforeAfterDrop_perm(:,6);
BeforeAfterDrop_perm(:,12)=BeforeAfterDrop_perm(:,10)+BeforeAfterDrop_perm(:,9);

figure; set(gcf,'position',[440   378   303   420]);
hold on;
Pos=1; data=BeforeAfterDrop(:,12); widthLine=2; widthBar=1.6; sizeDot=400; markerType='o';
colorBar=[.8 .8 .8;0 0 0];
line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)

patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
    [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');

xspread=(rand(1,length(data))-0.5)*widthBar/8+Pos-0.4;
scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);
 
            
Pos=2; data=BeforeAfterDrop_perm(:,12); widthLine=2; widthBar=1.6; sizeDot=400; markerType='o';
colorBar=[.4 .4 .4;0 0 0];
line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,25),'Color',colorBar(2,:),'LineWidth',widthLine)
line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*prctile(data,75),'Color',colorBar(2,:),'LineWidth',widthLine)
line([Pos-0.1*widthBar Pos-0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)
line([Pos-0.1*widthBar Pos+0.1*widthBar],[1 1].*nanmean(data),'Color',colorBar(2,:),'LineWidth',widthLine+1)
line([Pos+0.1*widthBar Pos+0.1*widthBar],[prctile(data,25) prctile(data,75)],'Color',colorBar(2,:),'LineWidth',widthLine)

patch([Pos-0.1*widthBar Pos+0.1*widthBar Pos+0.1*widthBar Pos-0.1*widthBar Pos-0.1*widthBar]',...
    [prctile(data,25) prctile(data,25) prctile(data,75) prctile(data,75) prctile(data,25)]',colorBar(1,:),'FaceAlpha',0.5,'EdgeColor','none');

xspread=(rand(1,length(data))-0.5)*widthBar/8+Pos-0.4;
scatter(xspread,data,'Marker',markerType,'MarkerFaceColor',colorBar(1,:),'MarkerEdgeColor',colorBar(2,:),'MarkerFaceAlpha',0.5,'SizeData',sizeDot/4);

format_fig;
set(gca,'XTick',1:2,'XTickLabel',{'real','perm'});
xlim([0.2 2.5])
ylabel('Proportion Sleep After Drop')

export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_DeltaPower_StageAfterDrop.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_DeltaPower_StageAfterDrop.png'],'-r 300')
