%-----------------------------------------------------------------
% Analyses EDISON - EEG Data
%------------------------------------------------------------------

clear all; close all; clc;

%%%%%% Init %%%%%%%
pathT = 'C:\Users\OUDIETTE.Delphine\Desktop\WagnerEdison project\Analyses_Results\T\';
SavePath = 'C:\Users\OUDIETTE.Delphine\Desktop\WagnerEdison project\Analyses_Results\Edison_Tables\';

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
    
    LastEpoch = T.Epoch(EndEdison)-1;
    if SleepEdison ==1
        LastStage = unique(T.Stade(T.Epoch ==LastEpoch));
        place_sleep = find(T.Stade>0);
        LastStage_Sleep =  T.Stade(place_sleep(end))
        Time_Wake_aftSleep = length(place_sleep(end):EndEdison);
    else
        LastStage = NaN;
        LastStage_Sleep = NaN;
        Time_Wake_aftSleep = NaN;
    end
    
    % Number of stage transition
    transition = diff(HypnoEdison);
    Nb_Stage_Transition = nnz(transition);
    
    % Micro-sleep transition
    N1_period = T.N1(StartEdison:EndEdison)';
    micro_transition = diff(N1_period);
    Nb_Micro_Sleep_Transition = nnz(micro_transition);
    
    % N1 & N2 epochs
    N1Epochs = T.Epoch(cellstr(T.StadeL)=="1");
    N1EpochsEdison = N1Epochs(N1Epochs>= StartEdisonEpoch & N1Epochs<=EndEdisonEpoch);
    N1EpochsEdison = unique(N1EpochsEdison); NbEpochsN1Edison = length(N1EpochsEdison);
    BERN_Epochs = T.Epoch(T.N1(:,1)==1);
    BERN_EpochsEdison = BERN_Epochs(BERN_Epochs>= StartEdisonEpoch & BERN_Epochs<=EndEdisonEpoch);
    BER_EpochsEdison = unique(BERN_EpochsEdison); NbEpochsBERN_Edison = length(BERN_EpochsEdison);
    if NbEpochsBERN_Edison ==0
        BERN_latency = NaN;
    else
        BERN_latency = minus(BER_EpochsEdison(1), StartEdisonEpoch)*30; % in secs % at least two epochs in a row
    end
    
    N2Epochs = T.Epoch(cellstr(T.StadeL)=="2");
    N2EpochsEdison = N2Epochs(N2Epochs>= StartEdisonEpoch & N2Epochs<=EndEdisonEpoch);
    N2EpochsEdison = unique(N2EpochsEdison); NbEpochsN2Edison = length(N2EpochsEdison);
    
    WakeEpochs = T.Epoch(cellstr(T.StadeL)=="W");
    WakeEpochsEdison = WakeEpochs(WakeEpochs>= StartEdisonEpoch & WakeEpochs<=EndEdisonEpoch);
    WakeEpochsEdison = unique(WakeEpochsEdison); NbEpochsWakeEdison = length(WakeEpochsEdison);
    
    if NbEpochsN1Edison ==0 && NbEpochsN2Edison ==0
        epochEdisonlatency =NaN;
        Edisonlatency = NaN;
    else
        if isempty(N1EpochsEdison) ==1 %take into account sub directly falling into N2
            epochEdisonlatency = N2EpochsEdison(1);
        else
            epochEdisonlatency = N1EpochsEdison(1);
        end
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
        timeN2bfBall = NaN;
        timeN1bfFirstBall = NaN;
        timeBallfromStart = NaN;
        lengthLastBout = NaN;
        NbboutsBeforeBall = NaN;
        Ball2End = NaN;
        timeFirstBallfromStart = NaN;
        NbboutsBeforeFirstBall = NaN;
        timeAASMN1bfBall = NaN;
    else
        Fall = nnz(T.Ball(:,1));
        timeBall=find(T.Ball==1);
        %time from start to ball drop
        
        for i = 1:length(timeBall)
            timeFirstBallfromStart = timeBall(1) - StartEdison;
            timeBallfromStart(i) = timeBall(i) - StartEdison;
            if isempty(StartBoutsN1Edison) ==1 % pas de périodes de N1
                lengthLastBout(i) = 0;   NbboutsBeforeBall(i) =0;
                NbboutsBeforeFirstBall= 0;
                
            elseif minus(timeBallfromStart(i), StartBoutsN1Edison) <0 % N1 après lâcher bouteille
                lengthLastBout(i)=0;NbboutsBeforeBall(i)=0;
                NbboutsBeforeFirstBall= 0;
            else
                [First(1) NbboutsBeforeFirstBall]= min(abs(timeFirstBallfromStart-StartBoutsN1Edison));
                [Index(i) NbboutsBeforeBall(i)]= min(abs(timeBallfromStart(i)-StartBoutsN1Edison));
                lengthLastBout(i) = LengthBoutsN1Edison(NbboutsBeforeBall(i));
                if i>1; NbboutsBeforeBall(i) = NbboutsBeforeBall(i) - NbboutsBeforeBall(i-1) ; end
                if lengthLastBout(i)<0; lengthLastBout(i)=0;NbboutsBeforeBall(i)=0;end
            end
        end
        
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
            timeAASMN1bfBall = nnz(T.Stade(StartEdison:timeBall));
            timeN2bfBall = nnz(T.Stade(StartEdison:timeBall,1)==2);
            timeN1bfFirstBall = timeN1bfBall;
        else
            timeBallFell(1) = timeBall(1) - StartEdison; %time from start of the first drop
            timeN1bfBall(1) = nnz(T.N1(StartEdison:timeBall(1),1)); % 1er fall: time of N1 from the beginning
            timeN2bfBall(1) = nnz(T.Stade(StartEdison:timeBall,1)==2);
            timeN1bfFirstBall = nnz(T.N1(StartEdison:timeBall(1),1));
            for i = 2:Fall %other fall time of N1 from the preceding fall
                timeBallFell(i) = timeBall(i) - timeBall(i-1); %time from start of the first drop
                timeN1bfBall(i) = nnz(T.N1(timeBall(1):timeBall(i),1));
                timeN2bfBall(i) = nnz(T.Stade(timeBall(1):timeBall(i),1)==2);
                timeAASMN1bfBall(i) = nnz(T.Stade(timeBall(1):timeBall(i)));
            end
        end
        
    end
    
    timeBallFell = {timeBallFell}; stadeBall = {stadeBall};timeN1bfBall = {timeN1bfBall}; timeN2bfBall = {timeN2bfBall};
    LengthBoutsN1Edison = {LengthBoutsN1Edison}; timeAASMN1bfBall = {timeAASMN1bfBall};
    lengthLastBout = {lengthLastBout}; NbboutsBeforeBall = {NbboutsBeforeBall};timeBallfromStart = {timeBallfromStart};
    Sujet = cellstr(sujet);
    
    Edison_EEG(i_sub,:) = table(Sujet,SleepEdison,HypnoEdisonWake,HypnoEdisonN1, HypnoEdisonN2, NbEpochsN2Edison,timeFirstBallfromStart, NbEpochsN1Edison,NbEpochsWakeEdison, BERN_latency, Edisonlatency, ...
        PercentN1Edison_Micro, timeN1Edison, Fall,timeAASMN1bfBall, timeBallFell, stadeBall, timeN1bfBall, timeN2bfBall,timeN1bfFirstBall, timeBallfromStart, lengthLastBout,NbboutsBeforeFirstBall, NbboutsBeforeBall,...
        TimeBtwLastN1_2End,NbWakefromN1Bouts, Nb_Stage_Transition,Nb_Micro_Sleep_Transition, LengthBoutsN1Edison, LastStage, LastStage_Sleep, Time_Wake_aftSleep);
    
    %         save(strcat(pathT,'Details_Edison', sujet), 'LengthBoutsN1Edison');
    
    keep('Edison_EEG', 'SavePath', 'pathT', 'AllSubjects');
end
save(strcat(SavePath,'Edison_EEG'), 'Edison_EEG');

% Lauch Behaviour Analyses if desired

%Analyses_Edison_Behaviour_Clean;
