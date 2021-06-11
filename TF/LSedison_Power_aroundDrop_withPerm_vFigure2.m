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

data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).

%% Extract all drop timings
all_time_drops=[];
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
BeforeDrop_perm=[];
AfterDrop_perm=[];

BeforeDrop=[];
AfterDrop=[];

StageDrop=[];

new_table=[];

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
        hdr=ft_read_header([Folder_Name filesep File_Name]);
        events=ft_read_event([Folder_Name filesep File_Name]);
        dat=ft_read_data([Folder_Name filesep File_Name]);
        
        drop_indexes=T.Epoch(find(T.Ball(:,1)))-1;
        
        for k=1 %:length(drops)
            this_line=match_str(data_clean.ID,SubdID);
            if isempty(this_line)
                continue;
            end
            time_drop=T.TimeID(drops(k));
            this_stage_drop=max(ASSM_Score(find(ASSM_Score_Time>=time_drop-30 & ASSM_Score_Time<=time_drop))); %max(unique(T.Stade(T.Epoch==drop_indexes(k) | T.Epoch==drop_indexes(k)+1)));
            
%             figure;
%             plot(ASSM_Score_Time,ASSM_Score)
%             hold on
%             scatter(time_drop,ASSM_Score(ASSM_Score_Time==time_drop))
%             title(SubdID)
%             pause;
            time_eeg=1/hdr.Fs:1/hdr.Fs:hdr.nSamples/hdr.Fs;
            if (time_drop-60)<Beg_Task
                continue;
            end
            %             temp_eeg=dat(3,time_eeg>=Beg_Task & time_eeg<=End_Task);
            temp_eeg_before=dat(3,time_eeg>=time_drop-50 & time_eeg<time_drop-2);
            temp_eeg_after=dat(3,time_eeg>time_drop+2 & time_eeg<time_drop+50);
            w_window=6*hdr.Fs;
            if size(temp_eeg_before,2)<w_window || size(temp_eeg_after,2)<w_window
                continue;
            end
            
            w_overlap=w_window/2;
            df=0.2;
            freqV=1:0.2:30;
            [pow_before,faxis] = pwelch(temp_eeg_before,w_window,w_overlap,freqV,hdr.Fs,'psd');
            [pow_after,faxis] = pwelch(temp_eeg_after,w_window,w_overlap,freqV,hdr.Fs,'psd');
            BeforeDrop=[BeforeDrop ; pow_before];
            AfterDrop=[AfterDrop ; pow_after];
            StageDrop=[StageDrop ; this_stage_drop];
            new_table=[new_table ; data_clean(this_line,:)];
            
            nPerm=0;
            fprintf('%2.0f%%\n',0)
            temp_pow=[];
            temp_pow2=[];
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
                
                time_eeg=1/hdr.Fs:1/hdr.Fs:hdr.nSamples/hdr.Fs;
                if (this_rand_time-60)<Beg_Task
                    continue;
                end
                %             temp_eeg=dat(3,time_eeg>=Beg_Task & time_eeg<=End_Task);
                temp_eeg_before=dat(3,time_eeg>this_rand_time-60 & time_eeg<this_rand_time);
                temp_eeg_after=dat(3,time_eeg>this_rand_time & time_eeg<=End_Task);
                
                if size(temp_eeg_before,2)<w_window || size(temp_eeg_after,2)<w_window
                    continue;
                end
                w_window=6*hdr.Fs;
                w_overlap=w_window/2;
                df=0.2;
                freqV=1:0.2:30;
                [pow_before,faxis] = pwelch(temp_eeg_before,w_window,w_overlap,freqV,hdr.Fs,'psd');
                [pow_after,faxis] = pwelch(temp_eeg_after,w_window,w_overlap,freqV,hdr.Fs,'psd');
                temp_pow=[temp_pow ; pow_before];
                temp_pow2=[temp_pow2 ; pow_after];
                
            end
            BeforeDrop_perm=[BeforeDrop_perm ; mean(temp_pow)];
            AfterDrop_perm=[AfterDrop_perm ; mean(temp_pow2)];
            
        end
    end
end


%% Delta power before probes
% figure
% plot(faxis, 
%%
figure;
plot(faxis, mean(10*log10(BeforeDrop)),'r')
hold on
plot(faxis, mean(10*log10(BeforeDrop_perm)),'k')
% plot(faxis, mean(log10(BeforeDrop_perm)),'r--')
% plot(faxis, mean(log10(AfterDrop_perm)),'b--')

%%
ColorsGroup=[178,171,210;
    253,184,99;
    230,97,1]/256;

DeltaBefore=mean(log10(BeforeDrop(:,faxis>1 & faxis<4)),2);
DeltaAfter=mean(log10(AfterDrop(:,faxis>1 & faxis<4)),2);
% 
DeltaBefore_perm=mean(log10(BeforeDrop_perm(:,faxis>1 & faxis<4)),2);
DeltaAfter_perm=mean(log10(AfterDrop_perm(:,faxis>1 & faxis<4)),2);

figure; set(gcf,'position',[440   378   303   420]);
hold on;
Pos=1; data=DeltaBefore; widthLine=2; widthBar=1.6; sizeDot=400; markerType='o';
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
 
            
Pos=2; data=DeltaBefore_perm; widthLine=2; widthBar=1.6; sizeDot=400; markerType='o';
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
ylabel('Delta Power Before Drop')

export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_DeltaPower_AvBeforeDrop.eps'],'-r 300')
export_fig([pwd filesep '..' filesep 'FigMat' filesep 'Edison_TF_DeltaPower_AvBeforeDrop.png'],'-r 300')

%%
new_table.SleepGroup=nan(size(new_table,1),1);
new_table.SleepGroup(new_table.HypnoEdisonN1==0  & new_table.HypnoEdisonN2==0)=0;
new_table.SleepGroup(new_table.HypnoEdisonN1~=0  & new_table.HypnoEdisonN2==0)=1;
new_table.SleepGroup(new_table.HypnoEdisonN1~=0  & new_table.HypnoEdisonN2~=0)=2;
new_table.SleepGroup=categorical(new_table.SleepGroup);

new_table.PowDelta=squeeze(mean(log10(BeforeDrop(:,faxis>1 & faxis<4)),2));
new_table.PowTheta=squeeze(mean(log10(BeforeDrop(:,faxis>4 & faxis<8)),2));
new_table.PowAlpha=squeeze(mean(log10(BeforeDrop(:,faxis>8 & faxis<12)),2));

new_table.DropGroup=StageDrop;
new_table.DropGroup=categorical(new_table.DropGroup);

%%
mdl = fitglm(new_table,'InsightPost~1+DropGroup')

%%
figure;
plot(faxis,mean(zscore(log10(BeforeDrop(new_table.DropGroup=="0",:)),[],2)));
hold on
plot(faxis,mean(zscore(log10(BeforeDrop(new_table.DropGroup=="1",:)),[],2)));
plot(faxis,mean(zscore(log10(BeforeDrop(new_table.DropGroup=="2",:)),[],2)));

%%
ColorsGroup=[178,171,210;
    253,184,99;
    230,97,1]/256;
figure; hold on
bar(1,mean(new_table.InsightPost(new_table.DropGroup=="0")),'FaceColor',ColorsGroup(1,:),'EdgeColor',ColorsGroup(1,:),'LineWidth',3);
bar(2,mean(new_table.InsightPost(new_table.DropGroup=="1")),'FaceColor',ColorsGroup(2,:),'EdgeColor',ColorsGroup(2,:),'LineWidth',3);
bar(3,mean(new_table.InsightPost(new_table.DropGroup=="2")),'FaceColor',ColorsGroup(3,:),'EdgeColor',ColorsGroup(3,:),'LineWidth',3);


%%
figure;
subplot(1,2,1);
format_fig;
tempX=nanzscore(new_table.PowDelta);
tempY=new_table.InsightPost;

bins=prctile(tempX,0:33:100);
bin_values(1,1)=mean(tempX(tempX<bins(2)));
bin_values(1,2)=mean(tempY(tempX<bins(2)));
bin_values(2,1)=mean(tempX(tempX>=bins(2) & tempX<=bins(3)));
bin_values(2,2)=mean(tempY(tempX>=bins(2) & tempX<=bins(3)));
bin_values(3,1)=mean(tempX(tempX>bins(3)));
bin_values(3,2)=mean(tempY(tempX>bins(3)));

bin_values_sem(1,1)=sem(tempX(tempX<bins(2)));
bin_values_sem(1,2)=sem(tempY(tempX<bins(2)));
bin_values_sem(2,1)=sem(tempX(tempX>=bins(2) & tempX<=bins(3)));
bin_values_sem(2,2)=sem(tempY(tempX>=bins(2) & tempX<=bins(3)));
bin_values_sem(3,1)=sem(tempX(tempX>bins(3)));
bin_values_sem(3,2)=sem(tempY(tempX>bins(3)));

errorbar(bin_values(:,1),bin_values(:,2),-bin_values_sem(:,2)/2,+bin_values_sem(:,2)/2,-bin_values_sem(:,1)/2,+bin_values_sem(:,1)/2);

subplot(1,2,2);
format_fig;
tempX=nanzscore(new_table.PowTheta);
tempY=new_table.InsightPost;

bins=prctile(tempX,0:33:100);
bin_values(1,1)=mean(tempX(tempX<bins(2)));
bin_values(1,2)=mean(tempY(tempX<bins(2)));
bin_values(2,1)=mean(tempX(tempX>=bins(2) & tempX<=bins(3)));
bin_values(2,2)=mean(tempY(tempX>=bins(2) & tempX<=bins(3)));
bin_values(3,1)=mean(tempX(tempX>bins(3)));
bin_values(3,2)=mean(tempY(tempX>bins(3)));

bin_values_sem(1,1)=sem(tempX(tempX<bins(2)));
bin_values_sem(1,2)=sem(tempY(tempX<bins(2)));
bin_values_sem(2,1)=sem(tempX(tempX>=bins(2) & tempX<=bins(3)));
bin_values_sem(2,2)=sem(tempY(tempX>=bins(2) & tempX<=bins(3)));
bin_values_sem(3,1)=sem(tempX(tempX>bins(3)));
bin_values_sem(3,2)=sem(tempY(tempX>bins(3)));

errorbar(bin_values(:,1),bin_values(:,2),-bin_values_sem(:,2)/2,+bin_values_sem(:,2)/2,-bin_values_sem(:,1)/2,+bin_values_sem(:,1)/2);
