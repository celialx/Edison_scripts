%-----------------------------------------------------------------
% Analyses EDISON - EEG Data
%------------------------------------------------------------------

clear all; close all; clc;

%%%%%% Init %%%%%%%
pathT = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\T\';
SavePath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Edison_Tables\';

% Rajouter 13AB & 66LY for data bottle (did not finish post but did breaks)

AllSubjects = {'01HC'; '02ZB'; '03JF'; '04BS'; '05MB'; '06NJ'; '07LB'; '08CB'; '09FN'; '10OL'; '11BD'; '12TG'; '14CD'...
    ; '15CV' ; '16LK'; '17ML'; '18LW'; '19HJ'; '20JJ'; '21IH'; '22HK';'24MK'; '25DS'; '26DB'; '27JS'; '28MD'; '29AN'...
    ; '30LM';'32DP'; '33NM'; '34AP'; '35NM'; '36DF'; '37SA'; '39NL'; '40EE'; '41RL'; '42AP'...
    ;'43AC'; '44LM'; '45AK'; '46SR'; '47OK'; '48KD'; '49LL'; '50CV'; '51EG'; '52NF'; '53EI'; ...
    '55BB' ; '56AB'; '57LL'; '58BP'; '59SC'; '60HD'; '61JF'; '62IB'; '63EG'; '64LC';...
    '67VF'; '68CY'; '69AM'; '70MA'; '71KN'; '72IG'; '73GF'; '74LW';'75CC'; '76LC';...
    '77JD'; '78BG'; '79CR'; '80MB';'81LB'; '82RD'; '83MB'; '84RD'; '85AA'; ...
    '87MD';'88SN'; '89AA'; '90ED'; '91MC'; '92GG'; '93RS'; '94CC'; '95FD'; '96CH';...
    '97AJ'; '98JB'; '99VB'; '100CP'; '101MS'; '102OB'; '103JZ'; '104MR'; '105OG';...
    '107AS'; '108JC'; '109JF'; '110SD'; '111MB'; '112GC'};

for i_sub=1:length(AllSubjects)
    
    sujet=AllSubjects{i_sub,1};
    
    FileName=strcat(pathT, filesep,'T_',sujet, '.mat');
    load(FileName);
    
    % Length of the break, first column, second is the duration of event
    StartEdison = find(T.Start(:,1)==1) ;StartEdisonEpoch= T.Epoch(StartEdison);
    EndEdison = find(T.End(:,1)==1); EndEdisonEpoch = T.Epoch(EndEdison);
    LengthEdison = length(StartEdison:EndEdison); % time break in secs
    NbEpochsEdison = length(unique(T.Epoch(StartEdison):T.Epoch(EndEdison)));
    
    % Info on hypnogram
    
    HypnoEdison = T.Stade(StartEdison:EndEdison);
    nantest= find(isnan(HypnoEdison)); HypnoEdison(nantest) =0; % sometimes, the end is missing
    HypnoEdisonWake = (nnz(~HypnoEdison)/length(HypnoEdison))*100; % Pourcentage Wake
    HypnoEdisonN1 = (length(find(HypnoEdison==1))/length(HypnoEdison))*100;
    HypnoEdisonN2 = (length(find(HypnoEdison==2))/length(HypnoEdison))*100;
    if HypnoEdisonN1 || HypnoEdisonN2 > 0 ; SleepEdison = 1; else SleepEdison = 0; end
    
    % N1 & N2 epochs
    N1Epochs = T.Epoch(cellstr(T.StadeL)=="1");
    N1EpochsEdison = N1Epochs(N1Epochs>= StartEdisonEpoch & N1Epochs<=EndEdisonEpoch);
    N1EpochsEdison = unique(N1EpochsEdison); NbEpochsN1Edison = length(N1EpochsEdison);
    
    N2Epochs = T.Epoch(cellstr(T.StadeL)=="2");
    N2EpochsEdison = N2Epochs(N2Epochs>= StartEdisonEpoch & N2Epochs<=EndEdisonEpoch);
    N2EpochsEdison = unique(N2EpochsEdison); NbEpochsN2Edison = length(N2EpochsEdison);
    
    
    if isempty(N1EpochsEdison) ==1
        epochEdisonlatency =NaN;
        Edisonlatency = NaN;
    else
        epochEdisonlatency = N1EpochsEdison(1);
        Edisonlatency = minus(epochEdisonlatency, StartEdisonEpoch)*30; % in secs % at least two epochs in a row
        % Edisonlatency = time it took from the beginning to get an epoch (30 sec)
        % of N1 in secs
    end
    
    %%%% Info on N1 - BERN scoring %%%%
    
    allN1Edison = sprintf('%d',[T.N1(StartEdison:EndEdison,1)]);
    EndBoutsN1Edison = strfind(allN1Edison, '10');
    NbWakefromN1Bouts = length(EndBoutsN1Edison);
    if isempty(EndBoutsN1Edison) ==1
        LastBoutsN1Edison = NaN;
        TimeBtwLastN1_2End = NaN;
    else
        LastBoutsN1Edison = EndBoutsN1Edison(end);
        TimeBtwLastN1_2End = minus(EndEdison,LastBoutsN1Edison);
    end
    StartBoutsN1Edison = strfind(allN1Edison, '01');
    NbBoutsN1Edison = numel(StartBoutsN1Edison); %nb bouts N1 Edison
    
    if isempty (EndBoutsN1Edison) ==1
        LengthBoutsN1Edison =NaN;
    else
        for b = 1:length(EndBoutsN1Edison)
            LengthBoutsN1Edison(b) = minus(EndBoutsN1Edison(b),  StartBoutsN1Edison(b));
        end
        
    end
    
    N1 = find(T.N1(:,1)==1);
    N1Edison = N1(N1>= StartEdison & N1<=EndEdison);
    timeN1Edison = length(N1Edison);  %in secs
    PercentN1Edison_Micro = (timeN1Edison/LengthEdison)*100; % Percent N1 with micro scoring
    
    %%%% Info on Ball %%%%
    
    % ball fall everytime subject in N1?
    
    if nnz(T.Ball(:,1)) ==0
        Fall = 0;
        timeBallFell= NaN;
        stadeBall = NaN;
        timeN1bfBall = NaN;
        timeBallfromStart = NaN;
        lengthboutBall = NaN;
        NbboutsBeforeBall = NaN;
        Ball2End = NaN;
    else
        Fall = nnz(T.Ball(:,1));
        timeBall=find(T.Ball==1);
        %time from start to ball drop
        
        
        timeBallfromStart = timeBall(1) - StartEdison;
        
        stadeN1Ball = T.N1(timeBall-2);
        AllstadeBall = T.Stade(timeBall-2);
        
        for b = 1:length(stadeN1Ball);
            if stadeN1Ball(b) ==1;
                stadeBall(b) =1;
            else
                stadeBall(b) = AllstadeBall(b);
            end
        end
        
        
        if Fall ==1
            timeBallFell = timeBall - StartEdison;
            timeN1bfBall = nnz(T.N1(StartEdison:timeBall,1));
        else
            timeBallFell(1) = timeBall(1) - StartEdison; %time from start of the first drop
            timeN1bfBall(1) = nnz(T.N1(StartEdison:timeBall(1),1)); % 1er fall: time of N1 from the beginning
            
            for i = 2:Fall %other fall time of N1 from the preceding fall
                timeBallFell(i) = timeBall(i) - timeBall(i-1); %time from start of the first drop
                timeN1bfBall(i) = nnz(T.N1(timeBall(1):timeBall(i),1));
            end
        end
        
        if stadeBall == 0
            lengthboutBall = NaN;
            NbboutsBeforeBall = NaN;
            Ball2End = minus(EndEdison, timeBall(end));
            
        else
            for u = 1:length (StartBoutsN1Edison)
                
                if Fall ==1
                    test(u) = minus(timeBallFell, StartBoutsN1Edison(u));
                elseif strcmp(sujet, '40EE')
                    test(u) = minus(timeBall(2)-StartEdison, StartBoutsN1Edison(u)); %ce sujet a des périodes de N1 après premier lacher de bouteille
                else
                    test(u) = minus(timeBallfromStart, StartBoutsN1Edison(u)); %réflechir à comment prendre en compte les autres drops ici
                end
            end
            
            Ball2End = minus(EndEdison, timeBall(end));
            [minValue, IndexBoutBall] = min(test(test>0));
            lengthboutBall = minus(timeBallfromStart, StartBoutsN1Edison(IndexBoutBall)); % nombre de secondes du bout de N1 précédant le lacher de balle
            NbboutsBeforeBall = IndexBoutBall-1; % nb bouts de N1 before lacher balle
            
        end
    end
    
    timeBallFell = {timeBallFell}; stadeBall = {stadeBall};timeN1bfBall = {timeN1bfBall}; LengthBoutsN1Edison = {LengthBoutsN1Edison};
    Sujet = cellstr(sujet);
    
    Edison_EEG(i_sub,:) = table(Sujet,SleepEdison,HypnoEdisonWake,HypnoEdisonN1, HypnoEdisonN2, NbEpochsN2Edison, NbEpochsN1Edison,Edisonlatency, ...
        PercentN1Edison_Micro, timeN1Edison, Fall, timeBallFell, stadeBall, timeN1bfBall, timeBallfromStart, lengthboutBall,NbboutsBeforeBall,...
        Ball2End,TimeBtwLastN1_2End, NbWakefromN1Bouts, LengthBoutsN1Edison);
    
    %         save(strcat(pathT,'Details_Edison', sujet), 'LengthBoutsN1Edison');
    
    keep('Edison_EEG', 'SavePath', 'pathT', 'AllSubjects');
    
end
save(strcat(SavePath,'Edison_EEG'), 'Edison_EEG');

% Lauch Behaviour Analyses if desired

Analyses_Edison_Behaviour;
