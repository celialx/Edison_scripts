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
for nF=1:length(files)
    File_Name=files(nF).name;
    Folder_Name=files(nF).folder;
    fprintf('... processing %s\n',File_Name);
    
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([Folder_Name filesep File_Name]);
    events=ft_read_event([Folder_Name filesep File_Name]);
    dat=ft_read_data([Folder_Name filesep File_Name]);
    
    %%% consider only the pause
    hypnogram = import_hypnogram([Folder_Name filesep File_Name(1:4) '_hypnogram.TXT']);
    scores=find(hypnogram.score~='?');
    start_pause=(scores(1)-1)*30*hdr.Fs+1;
    end_pause=scores(end)*30*hdr.Fs;
    
    temp_data=dat(1:3,(start_pause:end_pause));
    temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
    
    
    [twa_results]=twalldetectnew_TA_v2(temp_data,hdr.Fs,0);
    all_Waves=[];
    for nE=1:3
        all_Waves=[all_Waves ; [repmat([nF 1 nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
            cell2mat(twa_results.channels(nE).negzx)' ...
            cell2mat(twa_results.channels(nE).poszx)' ...
            cell2mat(twa_results.channels(nE).wvend)' ...
            cell2mat(twa_results.channels(nE).maxnegpk)' ...
            cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
            cell2mat(twa_results.channels(nE).maxpospk)' ...
            cell2mat(twa_results.channels(nE).maxpospkamp)' ...
            cell2mat(twa_results.channels(nE).mxdnslp)' ...
            cell2mat(twa_results.channels(nE).mxupslp)' ...
            cell2mat(twa_results.channels(nE).maxampwn)' ...
            cell2mat(twa_results.channels(nE).minampwn)' ...
            ]];
    end
    fprintf('\n')
    save([save_path filesep File_Name(1:end-4) '_allSW'],'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:3
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep File_Name(1:end-4) '_slowSW'],'slow_Waves','hdr','paramSW')
    
    %%% Compute ERPs
    if doERP
        for nEl=1:3
            temp_slow_Waves=slow_Waves(slow_Waves(:,3)==nEl,:);
            temp_ERP=[];
            for nW=1:size(temp_slow_Waves,1)
                start=temp_slow_Waves(nW,5);
                if start-0.5*hdr.Fs>1 && start+1.5*hdr.Fs<hdr.nSamples
                    temp_EEG=temp_data(nEl,(start-0.5*hdr.Fs):(start+1.5*hdr.Fs));
                    temp_EEG=temp_EEG-mean(temp_EEG(1:0.5*hdr.Fs));
                    temp_ERP=[temp_ERP ; temp_EEG];
                end
            end
            SW_ERP(nF,nEl,:)=mean(temp_ERP,1);
        end
    end
    
    %%% match with stages
    slow_Waves(:,16)=nan;
    MyStages={'W','1','2','3'};
    for nW=1:size(slow_Waves,1)
        onset=slow_Waves(nW,5);
        onset_segment=floor((slow_Waves(nW,5)+start_pause)/(30*hdr.Fs))+1;
        this_score=find(ismember(MyStages,hypnogram.score(onset_segment)));
        slow_Waves(nW,16)=this_score;
    end
    
    %%% compute densities
    for nEl=1:3
        for nStage=1:4
            SW_densities(nF,nEl,nStage)=sum(slow_Waves(:,3)==nEl & slow_Waves(:,16)==nStage)/sum(ismember(hypnogram.score,MyStages(nStage)))/2;
        end
        SW_amp(nF,nEl)=nanmean(slow_Waves(slow_Waves(:,3)==nEl,4));
        SW_freq(nF,nEl)=nanmean(1./((slow_Waves(slow_Waves(:,3)==nEl,7)-slow_Waves(slow_Waves(:,3)==nEl,5))/hdr.Fs));
    end
end

%%
figure;
hold on;
format_fig;
for nEl=1:3
    plot(-0.5:1/hdr.Fs:1.5,squeeze(mean(SW_ERP(:,nEl,:),1)),'LineWidth',2);
end
legend(hdr.label(1:3))
xlabel('Time from Wave onset (s)')
ylabel('Amplitude (\muV)')