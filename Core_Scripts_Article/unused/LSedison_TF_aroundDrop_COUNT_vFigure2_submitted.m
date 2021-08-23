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

ColorsGroup=[93 175 117;
    133 189 181;
    67 115 128]/256;

ColorsSolver=[196 146 46 ; 177 68 22]/256;
data_clean=readtable('../Edison_Tables/Clean_Data_SS.txt');

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).

%% Gather timing all drops
totPerm=100;
%%
all_Files=[];
doERP=1;
all_TF_drops=[];
StageDrop=[];
all_ERP=[];
all_maxAbs=[];
InsigthDrop=[];
nc=0;
totPerm=100;
all_TF_drops_perm=[];
StageDrop_perm=[];
PptionSleep_after=[];
PptionSleep_after_perm=[];
PptionSleep_before=[];
PptionSleep_before_perm=[];
all_keeps=[];
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
        threshArt           = 750;
        nc=nc+1;
        pption_rejTr(nc,:)  = [sum(maxAbs(:,3)>threshArt) size(maxAbs,1) 100*mean(maxAbs(:,3)>threshArt)];
        fprintf('... ... %g epochs (%g %%) rejected for threshold %g uV\n',sum(maxAbs(:,3)>threshArt),100*mean(maxAbs(:,3)>threshArt),threshArt)
        
        drops(maxAbs(:,3)>threshArt)=[];
        keep(maxAbs(:,3)>threshArt)=0;
        %          trl = LSedison_trialfun_drops(cfg)
        %         start_drop=(T.TimeID(find(T.Ball(:,1)))-1)*hdr.Fs;
                all_Files=[all_Files ; repmat({SubdID},1,1)];

        if isempty(find(keep))
         all_keeps=[all_keeps ; [nF 0 size(maxAbs,1) sum(maxAbs(:,3)>threshArt)]];
           continue
        else
                     all_keeps=[all_keeps ; [nF 1 size(maxAbs,1) sum(maxAbs(:,3)>threshArt)]];
        end
        
        
        
    end
end
%     load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')

%%
SumDrops=[];
for nF=1:length(all_Files)
    thisF=match_str(data_clean.Sujet,all_Files{nF});
    if ~isempty(thisF)
        if data_clean.InsightPre(thisF)
        SumDrop(nF)=NaN;
        else
        SumDrop(nF)=sum(~isnan(table2array(data_clean(thisF,51:54))));
        end
    else
        SumDrops(nF)=NaN;
    end
end
