%%
clear all
close all

path_raincloud='/Users/tand0009/WorkGit/projects/ext/RainCloudPlots/';
addpath(genpath(path_raincloud));

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

%% Extract all drop timings
all_time_drops=[];
totPerm=1000;
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
    %     figure;
    ASSM_Score=T.Stade;
    Beg_Task=T.TimeID(find(T.Start(:,1)));
    End_Task=T.TimeID(find(T.End(:,1)));
    ASSM_Score=ASSM_Score(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    ASSM_Score_Time=T.TimeID(T.TimeID>=Beg_Task & T.TimeID<=End_Task);
    
    drops=find(T.Ball(:,1));
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
plotdata=[];
for col = 1:2
    for plo = 1:2
        if plo==1
            plotdata{col, plo} = BeforeAfterDrop(:,4+3*col);
        elseif plo==2
            plotdata{col, plo} = BeforeAfterDrop_perm(:,4+3*col);
        end
    end
end

% make figure
cl=[];
cl(1,:)=[0 0 0.8];
cl(2,:)=[1 1 1]*0.7;

h   = rm_raincloud(plotdata, cl);

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

figure;
subplot(1,2,1)
format_fig;
hb=[];
simpleBarPlot(1+0.2,BeforeAfterDrop(:,6),[1 1 1;0 0 0.8],0.35,'r');
hb(1)=simpleBarPlot(2+0.2,BeforeAfterDrop(:,9),[0 0 0.8],0.35,'r');

simpleBarPlot(1-0.2,BeforeAfterDrop_perm(:,6),[1 1 1;0.8 0.8 0.8],0.35,'r');
hb(2)=simpleBarPlot(2-0.2,BeforeAfterDrop_perm(:,9),[0.8 0.8 0.8],0.35,'r');

xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Pre','Post'});
ylabel('Proportion of sleep')
legend(hb,{'Real','Perm'})
title('N1')

subplot(1,2,2)
format_fig;
hb=[];
simpleBarPlot(1+0.2,BeforeAfterDrop(:,7),[1 1 1;0 0 0.8],0.35,'r');
hb(1)=simpleBarPlot(2+0.2,BeforeAfterDrop(:,10),[0 0 0.8],0.35,'r');

simpleBarPlot(1-0.2,BeforeAfterDrop_perm(:,7),[1 1 1;0.8 0.8 0.8],0.35,'r');
hb(2)=simpleBarPlot(2-0.2,BeforeAfterDrop_perm(:,10),[0.8 0.8 0.8],0.35,'r');

xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Pre','Post'});
ylabel('Proportion of sleep')
legend(hb,{'Real','Perm'})
title('N2')

%%
BeforeAfterDrop2(:,11)=BeforeAfterDrop2(:,7)+BeforeAfterDrop2(:,6);
BeforeAfterDrop2(:,12)=BeforeAfterDrop2(:,10)+BeforeAfterDrop2(:,9);

BeforeAfterDrop_perm2(:,11)=BeforeAfterDrop_perm2(:,7)+BeforeAfterDrop_perm2(:,6);
BeforeAfterDrop_perm2(:,12)=BeforeAfterDrop_perm2(:,10)+BeforeAfterDrop_perm2(:,9);

figure;
subplot(1,2,1)
format_fig;
hb=[];
simpleBarPlot(1+0.2,BeforeAfterDrop2(:,6),[1 1 1;0 0 0.8],0.35,'r');
hb(1)=simpleBarPlot(2+0.2,BeforeAfterDrop2(:,9),[0 0 0.8],0.35,'r');

simpleBarPlot(1-0.2,BeforeAfterDrop_perm2(:,6),[1 1 1;0.8 0.8 0.8],0.35,'r');
hb(2)=simpleBarPlot(2-0.2,BeforeAfterDrop_perm2(:,9),[0.8 0.8 0.8],0.35,'r');

xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Pre','Post'});
ylabel('Proportion of sleep')
legend(hb,{'Real','Perm'})
title('N1')

subplot(1,2,2)
format_fig;
hb=[];
simpleBarPlot(1+0.2,BeforeAfterDrop2(:,7),[1 1 1;0 0 0.8],0.35,'r');
hb(1)=simpleBarPlot(2+0.2,BeforeAfterDrop2(:,10),[0 0 0.8],0.35,'r');

simpleBarPlot(1-0.2,BeforeAfterDrop_perm2(:,7),[1 1 1;0.8 0.8 0.8],0.35,'r');
hb(2)=simpleBarPlot(2-0.2,BeforeAfterDrop_perm2(:,10),[0.8 0.8 0.8],0.35,'r');

xlim([0.2 2.8])
set(gca,'XTick',1:2,'XTickLabel',{'Pre','Post'});
ylabel('Proportion of sleep')
legend(hb,{'Real','Perm'})
title('N2')