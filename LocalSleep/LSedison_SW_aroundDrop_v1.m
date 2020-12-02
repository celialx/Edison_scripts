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
all_PSTH_waves=cell(1,3);
StageDrop=[];
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
        
        load([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
        start_pause=(T.TimeID(find(T.Start(:,1)))-1)*hdr.Fs+1;
        start_drop=(T.TimeID(find(T.Ball(:,1)))-1)*hdr.Fs;
        stages_drop=T.Stade(find(T.Ball(:,1)));
        for k=1:length(start_drop)
            for nEl=1:3
                onset_waves=(slow_Waves(slow_Waves(:,3)==nEl,5)+start_pause-start_drop(k))/hdr.Fs;
                onset_waves=onset_waves(onset_waves>-2*60 & onset_waves<0.5*60);
                waves_bin=hist(onset_waves,-2*60:5:0.5*60);
                
                all_PSTH_waves{nEl}=[all_PSTH_waves{nEl} ; waves_bin];
            end
            StageDrop=[StageDrop stages_drop(k)];
        end
    end
end

%%
for nEl=1:3
    figure;
    times=-2*60:5:0.5*60;
    temp_toplot=all_PSTH_waves{nEl}(:,:);
    simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,0,1,1);
    temp_toplot=all_PSTH_waves{nEl}(StageDrop==0,:);
    simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,0,1,1);
    temp_toplot=all_PSTH_waves{nEl}(StageDrop~=0,:);
    simpleTplot(times,temp_toplot,0,'b',0,'-',0.5,1,0,1,1);
    xlim([-60 20])
    format_fig;
    xlabel('Time from Drop (s)')
    ylabel('Power (dB)')
    title(sprintf('N=%g drops',size(all_PSTH_waves{2},1)))
    line([0 0],ylim,'Color','k','LineStyle','--')
end
% %  hold on;
% %  temp_toplot=squeeze(mean(all_TF_drops(StageDrop==0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
% %  simpleTplot(times,temp_toplot,0,'r',0,'-',0.5,1,0,1,1);
% %  temp_toplot=squeeze(mean(all_TF_drops(StageDrop~=0,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
% %  simpleTplot(times,temp_toplot,0,'b',0,'-',0.5,1,0,1,1);
%