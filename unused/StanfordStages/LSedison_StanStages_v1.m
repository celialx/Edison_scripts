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

%% Table
load ../Edison_Tables/All_Data.mat

data_clean=Data(Data.InsightPre==0,:);
data_clean(:,[20 41 62 81 ])=[];

data_clean.tot_epochs=nan(size(data_clean,1),1);
data_clean.AASM_W=nan(size(data_clean,1),1);
data_clean.AASM_N1=nan(size(data_clean,1),1);
data_clean.AASM_N2=nan(size(data_clean,1),1);

data_clean.STAGES_W=nan(size(data_clean,1),1);
data_clean.STAGES_N1=nan(size(data_clean,1),1);
data_clean.STAGES_N2=nan(size(data_clean,1),1);
data_clean.STAGES_N3=nan(size(data_clean,1),1);
data_clean.STAGES_RE=nan(size(data_clean,1),1);

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
    
    this_line=match_str(data_clean.ID,SubdID);
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
    Beg_Task=T.Epoch(find(T.Start(:,1)));
    End_Task=T.Epoch(find(T.End(:,1)));
    ASSM_Score=ASSM_Score(T.Epoch>=Beg_Task & T.Epoch<=End_Task);
    ASSM_Score_Time=T.TimeID(T.Epoch>=Beg_Task & T.Epoch<=End_Task);
    
    
    hypnogram_perSec=reshape(repmat(hypnogram,1,15)',1,numel(repmat(hypnogram,1,15)));
    time_hypno_perSec=1:length(hypnogram_perSec);
    Stanford_Score=hypnogram_perSec(time_hypno_perSec>=(Beg_Task-1)*30+1 & time_hypno_perSec<=30*End_Task);
    Stanford_Score_Time=time_hypno_perSec(time_hypno_perSec>=(Beg_Task-1)*30+1 & time_hypno_perSec<=30*End_Task);
    Stanford_HypnoD=nan(5,size(Stanford_Score,2));
    for j=1:size(hypnodensity,2)
        temparray=table2array(hypnodensity(:,j));
        temparray2=reshape(repmat(temparray,1,15)',1,numel(repmat(temparray,1,15)));
        Stanford_HypnoD(j,:)=temparray2(time_hypno_perSec>=(Beg_Task-1)*30+1 & time_hypno_perSec<=30*End_Task);
    end
    
    nc=nc+1;
    perS_deepestSt(nc)=max(ASSM_Score);
    
    data_clean.tot_epochs(this_line)=length(ASSM_Score)/30;
    data_clean.AASM_W(this_line)=sum(ASSM_Score==0)/30;
    data_clean.AASM_N1(this_line)=sum(ASSM_Score==1)/30;
    data_clean.AASM_N2(this_line)=sum(ASSM_Score==2)/30;
    
    data_clean.STAGES_W(this_line)=nanmean(Stanford_HypnoD(1,:));
    data_clean.STAGES_N1(this_line)=nanmean(Stanford_HypnoD(2,:));
    data_clean.STAGES_N2(this_line)=nanmean(Stanford_HypnoD(3,:));
    data_clean.STAGES_N3(this_line)=nanmean(Stanford_HypnoD(4,:));
    data_clean.STAGES_RE(this_line)=nanmean(Stanford_HypnoD(5,:));
    
end

%%
writetable(data_clean,'../Edison_Tables/Clean_Data_SS.txt');

