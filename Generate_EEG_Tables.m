%% Pour ce script il faut être dans un Dossier Scripts et dans la racine de ce dossier avoir un dossier qui contient 
%les EDF, les stades et un dossier qui contient les Evenemnts Scorés. Il
%faut aussi avoir les fonctions suivantes dans le dossier script
%            - blockEdfLoad.m
%            - importStade.m
%            - importfileEvent.m
clear all


%------------------------------------------------------------------------------------------------------
% Path
%------------------------------------------------------------------------------------------------------

pathD = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\FinalScoringEdison\';
pathEEG = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\EEG\';
pathT = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\T\';

%------------------------------------------------------------------------------------------------------
% Initialisation
%------------------------------------------------------------------------------------------------------

% 13AB for supp data on the bottle (with subjects of Insight) but remove the subject for Edison analyses
% Il faut penser à convertir les fichiers evenementscores en csv!!

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

%Télécharge l'EDF
FileName=strcat(pathEEG,sujet,'_EEG.edf'); % Télécharge l'EDF
[header,signalHeader,signalCell] = blockEdfLoad(FileName);

%------------------------------------------------------------------------------------------------------
% Création d'une matrice T qui contient une ligne par seconde pour toute la
% sieste
%------------------------------------------------------------------------------------------------------

% Struture de la matrice à partir de la première seconde extraite de l'EDF
%----------------------------------------------------------------------------

% Transforme le temps initial en variable h:m:s pour pouvoir remplir la matrice
h=str2double(header.recording_starttime(1:2));
m=str2double(header.recording_starttime(4:5));
s=str2double(header.recording_starttime(7:8));

% Longueur de la matrice en seconde : 
lenT = header.num_data_records;

% Ici on définit toutes les variables que contriendra la matrice T
Time          = cell(lenT,1); % Temps en h:m:s
TimeID        = zeros(lenT,1); % Temps en n°
Epoch         = zeros(lenT,1); % Le N° de l'epoch de 30 sec
Stade         = zeros(lenT,1); % Le stade de la seconde en numero
StadeL        = cell(lenT,1); % Le stade de la seconde en lettre

N1    = zeros(lenT,1); % L'evenement 1 qui est exporté du fichier 
Start    = zeros(lenT,1); % L'evenement 2 qui est exporté du fichier 
End    = zeros(lenT,1); % L'evenement 3 qui est exporté du fichier 
Ball    = zeros(lenT,1); % L'evenement 4 qui est exporté du fichier 
% Probes    = zeros(lenT,1); % L'evenement 5 qui est exporté du fichier 

SR=signalHeader.samples_in_record; % Sample rate


% Cette boucle remplit les variables Time et Time ID en débutant depuis la
% première seconde de l'enregristrement extraite grâce à l'EDF
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

% Incorpore les stades à la matrice
%----------------------------------------------------------------------------

%Téléchargement du score pour l'incorporer à la matrice
filename = strcat(pathD,sujet,'_hypnogram.txt');
StadeImported = importStade(filename);
StadeImported{end+1,1}=StadeImported{end,1}; % Il manque souvent un stade à la fin on le rajoute et il est défini comme le même que celui d'avant, souvent de l'éveil

%Stade Importer est une cellule, les stades sont sous forme de texte, cette
%bloucle permet de les mettre sous la forme de nombre
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

%Remplit les variables Stades Et Stades L grace aux stades extrais comme 
%chaque époque fait 30 secondes, on rempli avec epoch 1 pour 30 secondes
%puis epoch 2 pour 30 secondes etc...
j=1;
    for i_T = 1:lenT
        Epoch(i_T,1)=j;
        Stade(i_T,1)=StadeNum(j,1);
        StadeL{i_T,1}=StadeImported{j,1};

        if mod(i_T,30)==0
            j=j+1;
        end
    end
    
% Incorpore les évenements à la matrice depuis le fichier évement
%----------------------------------------------------------------------------
filename = strcat(pathD,sujet,'_evenementsscores.csv');
ScoredEvents = importfileEvent(filename);

Text={'1';'2';'3';'4'};
nbrevent=length(Text); % Nombre d'evenements possibles

Mat01      = {N1;Start;End;Ball}; %Matice qui contient un 1 si il y a un evement à cette seconde

Duration   = zeros(lenT,nbrevent); % Duration of event ?

Evenements =struct('Text',Text,'Mat01',Mat01);

%Là il faut traduire les evenements en 0 et 1 dans chaques catégories 
for i_se=1:height(ScoredEvents)
    
    i_time=find(strcmp(ScoredEvents.Time{i_se,1},Time)); % Trouve à quelle seconde s'applique l'evenement

    for i_evt=1:nbrevent  % cherche quel evenement c'est
        if sum(ismember(ScoredEvents.Event{i_se,1},Evenements(i_evt).Text))==1 %Cherche  
            Evenements(i_evt).Mat01(i_time,1)=Evenements(i_evt).Mat01(i_time,1)+1;

            D=ScoredEvents.Duree{i_se,1};
            % Sachant que dans la matrice il y a des durée en seconde : 0.1
            % ou en millisecondes il faut prendre différement les deux cas
            % c'est ce que fait la commande en dessous
            if length(ScoredEvents.Duree{i_se,1})==5 
                Dnum=str2double(D(1,end-1:end));
            else
                Dnum=str2double(D(1,end-3:end));
            end


            if Dnum > 1.5 && i_evt ~= 2 && i_evt ~= 3  % Si un évenement dure plus qu'une seconde on le marque sur les secondes d'après
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