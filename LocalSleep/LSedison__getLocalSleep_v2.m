%%
clear all
close all

Computer = 'Célia';
if strcmp(Computer, 'Thomas')
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
path_localsleep='/Users/tand0009/WorkGit/projects/inprogress/wanderIM/localsleep';
addpath(path_fieldtrip);
path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
data_path='/Users/tand0009/Data/LS_Edison/';
save_path='/Users/tand0009/Data/LS_Edison/LocalSleep';
elseif strcmp(Computer, 'Célia')
    path_fieldtrip='D:\MATLAB\Toolbox\fieldtrip-20200409\';
        path_localsleep='D:\MATLAB\Toolbox\wanderIM\localsleep\';
    path_LSCPtools='D:\MATLAB\Toolbox\LSCPtools\';
save_path='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\LocalSleep\';
        data_path='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\EEG\';
T_path = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\T\';
path_signal = 'D:\MATLAB\Toolbox\signal\';
end

addpath(path_localsleep);
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath(path_signal));

files=dir([data_path '*.edf']);

%% INFO FROM CELIA
% - le fichier EDF: attention on enregistrait l'intégralité de l'expérience donc la pause ne commence normalement pas au début de l'enregistrement (sauf si on a oublié de lancer l'enregistrement au tout début...!)
% - un fichier mat avec une matrice seconde par seconde de l'enregistrement. Je pense que ce qui est important pour toi, c'est les colonnes Start et End. A chaque fois, le premier 1 de chacune des colonnes représente le début et la fin de la pause. Tu souhaites donc regarder les tracés EEG entre ces deux marqueurs.
% - Un fichier excel qui te redonne ces mêmes infos de manière peut-être plus rapide et visuelle (le moment du start, l'epoch, le scoring etc).


%%
doERP=1;
for nF=1:length(files)
    File_Name=files(nF).name;
%     Folder_Name=files(nF).folder;
     Folder_Name=data_path;
    SubdID=File_Name;
    sep=findstr(File_Name,'_');
    SubdID=SubdID(1:sep(1)-1);
    
    if exist([T_path 'T_' SubdID '.mat'])==0
        warning('matrix file missing');
        continue;
    end
    load([T_path 'T_' SubdID '.mat'])
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([Folder_Name filesep File_Name]);
    events=ft_read_event([Folder_Name filesep File_Name]);
    dat=ft_read_data([Folder_Name filesep File_Name]);
    
    %%% consider only the pause
    %     hypnogram = import_hypnogram([Folder_Name filesep File_Name(1:4) '_hypnogram.TXT']);
    %     scores=find(T.StadeL~='?');
    start_pause=(T.TimeID(find(T.Start(:,1)))-1)*hdr.Fs+1;
    end_pause=(T.TimeID(find(T.End(:,1))))*hdr.Fs;
    
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
                if start-0.5*hdr.Fs>1 && start+1.5*hdr.Fs<size(temp_data,2)
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
    MyStagesCode=[0 1 2 3];
    for nW=1:size(slow_Waves,1)
        onset=slow_Waves(nW,5);
        onset_segment=floor((slow_Waves(nW,5)+start_pause)/(hdr.Fs));
        this_score=find(ismember(MyStagesCode,T.Stade(T.TimeID==onset_segment)));
        if isempty(this_score)
        slow_Waves(nW,16)=nan;
        else
        slow_Waves(nW,16)=this_score;
        end
    end
    
    %%% compute densities
    for nEl=1:3
        for nStage=1:4
            SW_densities(nF,nEl,nStage)=sum(slow_Waves(:,3)==nEl & slow_Waves(:,16)==nStage)/sum(ismember(T.Stade,MyStagesCode(nStage)))/2;
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

%%
figure;
Colors=[0.4000    0.7608    0.6471;0.9882    0.5529    0.3843;0.5529    0.6275    0.7961];
for nEl=1:3
    for nSta=1:3
        temp=squeeze(SW_densities(:,nEl,nSta))*2*60;
        simpleBarPlot(nEl+(nSta-2)*.2,temp,[Colors(nSta,:)],0.18,'k');
    end
end
format_fig;
set(gca,'XTick',1:3,'XTickLabel',{'Fp1','C3','O1'})
ylabel('SW density (wave/min)')
xlim([0.2 3.8])