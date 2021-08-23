% Script to generate T matrix (one-by-one second table with sleep scoring
% and events like bottle drops)

% You need: EDF files, scoring files and evenementscores files, as well
% as the following scripts (present in the folder
% "Pack_Generate_EEG_Tables"):
%            - blockEdfLoad.m
%            - importStade.m
%            - importfileEvent.m

clear all; close all; clc;

%------------------------------------------------------------------------------------------------------
% Path
%------------------------------------------------------------------------------------------------------

pathD = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\FinalScoringEdison\';
pathEEG = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\EEG\';
pathT = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\T\';

%------------------------------------------------------------------------------------------------------
% Initialisation
%------------------------------------------------------------------------------------------------------

% The "evenementscores" files need to be converted into csv first.

AllSubjects = {'01HC'; '02ZB'; '03JF'; '04BS'; '05MB'; '06NJ'; '07LB'; '08CB'; '09FN'; '10OL'; '11BD'; '12TG'; '14CD'...
   ; '15CV' ; '16LK'; '17ML'; '18LW'; '19HJ'; '20JJ'; '21IH'; '22HK';'24MK'; '25DS'; '26DB'; '27JS'; '28MD'; '29AN'...
    ; '30LM'; '31MS'; '32DP'; '33NM'; '34AP'; '35NM'; '36DF'; '37SA'; '39NL'; '40EE'; '41RL'; '42AP'...
    ;'43AC'; '44LM'; '45AK'; '46SR'; '47OK'; '48KD'; '49LL'; '50CV'; '51EG'; '52NF'; '53EI'; ...
    '55BB' ; '56AB'; '57LL'; '58BP'; '59SC'; '60HD'; '61JF'; '62IB'; '63EG'; '64LC';...
    '67VF'; '68CY'; '69AM'; '70MA'; '71KN'; '72IG'; '73GF'; '74LW';'75CC'; '76LC';...
    '77JD'; '78BG'; '79CR'; '80MB';'81LB'; '82RD'; '83MB'; '84RD'; '85AA'; ...
    '87MD';'88SN'; '89AA'; '90ED'; '91MC'; '92GG'; '93RS'; '94CC'; '95FD'; '96CH';...
    '97AJ'; '98JB'; '99VB'; '100CP'; '101MS'; '102OB'; '103JZ'; '104MR'; '105OG';...
    '107AS'; '108JC'; '109JF'; '110SD'; '111MB'; '112GC'};

for i_sub=1:length(AllSubjects)

sujet=AllSubjects{i_sub,1};

% Load EDF
FileName=strcat(pathEEG,sujet,'_EEG.edf');
[header,signalHeader,signalCell] = blockEdfLoad(FileName);

%------------------------------------------------------------------------------------------------------
% Create a matrix T that contains one ligne per second for the entirety of
% the nap
%------------------------------------------------------------------------------------------------------

% Transform initial time into h:m:s variable to fill the matrix 
h=str2double(header.recording_starttime(1:2));
m=str2double(header.recording_starttime(4:5));
s=str2double(header.recording_starttime(7:8));

% Matrix length in second
lenT = header.num_data_records;

% Define all matrix variables
Time          = cell(lenT,1); % Time in h:m:s
TimeID        = zeros(lenT,1); % Time in n°
Epoch         = zeros(lenT,1); % 30s-epoch number
Stade         = zeros(lenT,1); % Sleep scoring
StadeL        = cell(lenT,1); % Sleep scoring (letter)

N1    = zeros(lenT,1); % Event 1 from scoring file
Start    = zeros(lenT,1); % Event 2 from scoring file
End    = zeros(lenT,1); % Event 3 from scoring file
Ball    = zeros(lenT,1); % Event 4 from scoring file

SR=signalHeader.samples_in_record; % Sample rate


% Loop filling variables Time and Time ID starting with first second
%extracted from the EDF file
for i = 1:header.num_data_records

    if s<10
        sstr=strcat('0',num2str(s));
    else
        sstr=num2str(s);
    end
    
    if m<10
        mstr=strcat('0',num2str(m));
    else
        mstr=num2str(m);
    end
    
    if h<10
        hstr=strcat('0',num2str(h));
    else
        hstr=num2str(h);
    end
        
    Time{i}=strcat(hstr,':',mstr,':',sstr);
    TimeID(i)=i;
    
    
    if s == 59 && m==59
        s=0;
        m=0;
        h=h+1;
    elseif s == 59
        s=0;
        m=m+1;
    else
       s=s+1; 
    end
   
end

% Integrate sleep substages to the matrix
%----------------------------------------------------------------------------

%Download score to integrate it into the matrix
filename = strcat(pathD,sujet,'_hypnogram.txt');
StadeImported = importStade(filename);
% The last stage is often missing, add the same than the previous one
StadeImported{end+1,1}=StadeImported{end,1};

%Loop to transform stade into number
nbrEpoch=height(StadeImported);

StadeNum=zeros(nbrEpoch,1);
for i_stadeimported = 1:nbrEpoch
        if strcmp(StadeImported{i_stadeimported,1},'W')==1
            StadeNum(i_stadeimported,1)=0;
        elseif strcmp(StadeImported{i_stadeimported,1},'R')==1
            StadeNum(i_stadeimported,1)=5;
        elseif strcmp(StadeImported{i_stadeimported,1},'?')==1
            StadeNum(i_stadeimported,1)=NaN;
        elseif strcmp(StadeImported{i_stadeimported,1},'1')==1
            StadeNum(i_stadeimported,1)=1;
        elseif strcmp(StadeImported{i_stadeimported,1},'2')==1
            StadeNum(i_stadeimported,1)=2;
        elseif strcmp(StadeImported{i_stadeimported,1},'3')==1
            StadeNum(i_stadeimported,1)=3;
        end
end

% Fill variables Stades and Stades L. One epoch = 30s
j=1;
    for i_T = 1:lenT
        Epoch(i_T,1)=j;
        Stade(i_T,1)=StadeNum(j,1);
        StadeL{i_T,1}=StadeImported{j,1};

        if mod(i_T,30)==0
            j=j+1;
        end
    end
    
% Integrate events into the matrix
%----------------------------------------------------------------------------
filename = strcat(pathD,sujet,'_evenementsscores.csv');
ScoredEvents = importfileEvent(filename);

Text={'1';'2';'3';'4'};
nbrevent=length(Text); % Nb of events

% Place a 1 at the corresponding second if there an event
Mat01      = {N1;Start;End;Ball}; 

Duration   = zeros(lenT,nbrevent); % Duration of event ?

Evenements =struct('Text',Text,'Mat01',Mat01);

% Place events into categories
for i_se=1:height(ScoredEvents)
    
    i_time=find(strcmp(ScoredEvents.Time{i_se,1},Time)); % Which second is the event

    for i_evt=1:nbrevent  % Search which event it is
        if sum(ismember(ScoredEvents.Event{i_se,1},Evenements(i_evt).Text))==1 
            Evenements(i_evt).Mat01(i_time,1)=Evenements(i_evt).Mat01(i_time,1)+1;

            D=ScoredEvents.Duree{i_se,1};
            % Consider whether the duration is in ms or s
            if length(ScoredEvents.Duree{i_se,1})==5 
                Dnum=str2double(D(1,end-1:end));
            else
                Dnum=str2double(D(1,end-3:end));
            end

% If an event lasts more than one second, put it on the following second
            if Dnum > 1.5 && i_evt ~= 2 && i_evt ~= 3
                for i_d = 1:round(Dnum)-1
                    Evenements(i_evt).Mat01(i_time+i_d,1)=1; 
                end
            end


                 Duration(i_time,i_evt)=Dnum;
        end
    end



    N1 = [Evenements(1).Mat01,Duration(:,1)];
    Start = [Evenements(2).Mat01,Duration(:,2)];
    End = [Evenements(3).Mat01,Duration(:,3)];
    Ball = [Evenements(4).Mat01,Duration(:,4)];

end

T=table(Time,TimeID,Epoch,Stade,StadeL,N1,Start,End,Ball);

save([pathT,'T_',sujet],'T');

end