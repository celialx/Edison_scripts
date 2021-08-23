% Comparison between DO and CL scorers

%------------------------------------------------------------------------------------------------------
% Path
%------------------------------------------------------------------------------------------------------
clear all;
close all;
clc;

savePath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\';
pathCL = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Two_scorers\T_CL\';
pathDO = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Two_scorers\T_DO\';
pathHypnoDO = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Two_scorers\DelphineScoringEdison\'
pathHypnoCL = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Two_scorers\CeliaScoringEdison\'

% Initialisation

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
    
    FileNameCL=strcat(pathCL, filesep,'T_CL_',sujet, '.mat'); % Télécharge l'EDF
    FileNameDO=strcat(pathDO, filesep,'T_DO_',sujet, '.mat'); % Télécharge l'EDF
    
    load(FileNameCL);  load(FileNameCL);
    
    % Length of the break
    StartEdison = find(T_CL.Start(:,1)==1); % Première colomne (la 2ème, c'est la durée de l'évènement)
    StartEdisonEpoch= T_CL.Epoch(StartEdison);
    
    EndEdison = find(T_CL.End(:,1)==1);
    EndEdisonEpoch = T_CL.Epoch(EndEdison);
    
    hypnoDO = strcat(pathHypnoDO,sujet,'_hypnogramDO.txt');
    hypnoCL = strcat(pathHypnoCL,sujet,'_hypnogramCL.txt');
    
    StadeImportedCL = table2array(importStade(hypnoCL)); StadeImportedCL= StadeImportedCL(StartEdisonEpoch:EndEdisonEpoch)
    StadeImportedDO= table2array(importStade(hypnoDO));StadeImportedDO= StadeImportedDO(StartEdisonEpoch:EndEdisonEpoch)
    
    
    for stade = 1:length(StadeImportedCL)
        if strcmp(StadeImportedCL{stade,1},'W')==1
            StadeImportedCL(stade,1)=0;
        elseif strcmp(StadeImportedCL{stade,1},'?')==1
            StadeImportedCL(stade,1)=NaN;
        elseif strcmp(StadeImportedCL{stade,1},'1')==1
            StadeImportedCL(stade,1)=1;
        elseif strcmp(StadeImportedCL{stade,1},'2')==1
            StadeImportedCL(stade,1)=2;
        end
    end
    
    
    
    for stade = 1:length(StadeImportedDO)
        if strcmp(StadeImportedDO{stade,1},'W')==1
            StadeImportedDO(stade,1)=0;
        elseif strcmp(StadeImportedDO{stade,1},'?')==1
            StadeImportedDO(stade,1)=NaN;
        elseif strcmp(StadeImportedDO{stade,1},'1')==1
            StadeImportedDO(stade,1)=1;
        elseif strcmp(StadeImportedDO{stade,1},'2')==1
            StadeImportedDO(stade,1)=2;
        end
    end
    
    
    
    
    Scoring(i_sub).CL = str2double(StadeImportedCL);
    Scoring(i_sub).DO = str2double(StadeImportedDO);
    
end

Scoring_T = struct2table(Scoring);
AllScoring= cell2mat(table2cell(Scoring_T));
xlswrite([savePath 'AllScoring.xlsx'], AllScoring);