%-----------------------------------------------------------------
% Analyses EDISON - Behaviour Data
%------------------------------------------------------------------

% The script "Analyses_Edison_EEG" need to be run before this one
% Think about adding mean PVT post for 84RD (bug script)

clear all; close all; clc;

% Paths
SavePath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Edison_Tables\';
DataPath = 'C:\Users\Célia\Desktop\WagnerEdison project\NRT_07_07_2020\Results\';
cd(DataPath);

AllSubjects = {'01HC'; '02ZB'; '03JF'; '04BS'; '05MB'; '06NJ'; '07LB'; '08CB'; '09FN'; '10OL'; '11BD'; '12TG'; '14CD'...
    ;'15CV' ; '16LK'; '17ML'; '18LW'; '19HJ'; '20JJ'; '21IH'; '22HK';'24MK'; '25DS'; '26DB'; '27JS'; '28MD'; '29AN'...
    ; '30LM';'32DP'; '33NM'; '34AP'; '35NM'; '36DF'; '37SA'; '39NL'; '40EE'; '41RL'; '42AP'...
    ;'43AC'; '44LM'; '45AK'; '46SR'; '47OK'; '48KD'; '49LL'; '50CV'; '51EG'; '52NF'; '53EI'; ...
    '55BB' ; '56AB'; '57LL'; '58BP'; '59SC'; '60HD'; '61JF'; '62IB'; '63EG'; '64LC';...
    '67VF'; '68CY'; '69AM'; '70MA'; '71KN'; '72IG'; '73GF'; '74LW';'75CC'; '76LC';...
    '77JD'; '78BG'; '79CR'; '80MB';'81LB'; '82RD'; '83MB'; '84RD'; '85AA'; ...
    '87MD';'88SN'; '89AA'; '90ED'; '91MC'; '92GG'; '93RS'; '94CC'; '95FD'; '96CH';...
    '97AJ'; '98JB'; '99VB'; '100CP'; '101MS'; '102OB'; '103JZ'; '104MR'; '105OG';...
    '107AS'; '108JC'; '109JF'; '110SD'; '111MB'; '112GC'};

% Load EEG info
EEGPath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Edison_Tables';
EEG_filename = [EEGPath filesep 'Edison_EEG.mat'];
load(EEG_filename);

% Load subject info (insight pre/post, hypnagogie etc)
SubInfoPath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results';
SubInfo_filename = [SubInfoPath filesep 'Subject_info.xlsx'];
Edison_Subject = readtable(SubInfo_filename);

% Check that we are comparing the same subject in each file
for n=1:length(AllSubjects)
    
    Subject_EEG = Edison_EEG.Sujet{n,1};
    Subject_Info = Edison_Subject.Code{n,1};
    Subject_Behav = AllSubjects{n,1};
    
    if ~strcmp(Subject_Behav, Subject_Info) || ~strcmp(Subject_Behav, Subject_EEG)
        error('not the same subject in info and matlab files');
    end
end

% Behaviour analyses
Nb_Pre_Trial = 60;
Nb_Post_Trial = 270; % without implicit

for i_sub=1:length(AllSubjects)
    
    sujet=AllSubjects{i_sub,1};
    Edison_Behaviour(i_sub).ID = sujet;
    
    Dir_Sub = get_subdir_regex(DataPath, sujet);
    Practice_Sub = get_subdir_regex_files(Dir_Sub, 'practice');
    Pre_Sub = get_subdir_regex_files(Dir_Sub, 'pre');
    Post_Sub = get_subdir_regex_files(Dir_Sub, 'post');
    
    % Check that we have all the data files (3 per subject)
    if strcmp(sujet, '25DS') % Pre solver, stopped after pre
        disp('hi');
    else
        if length(Practice_Sub) ~=1 || length(Pre_Sub) ~=1 || length(Post_Sub) ~=1
            error('Not the appropriate number of files');
        end
    end
    
    % Practice data
    load(cell2mat(Practice_Sub)); Practice_T = struct2table(Practice);
    if strcmp(sujet, '57LL') % Pb script, n'a enregistré que le dernier essai de practice
        Edison_Behaviour(i_sub).Correct_Practice =NaN; Edison_Behaviour(i_sub).RT_Practice = NaN;
        Edison_Behaviour(i_sub).TotalTime_Practice = NaN;
    else
        Edison_Behaviour(i_sub).Correct_Practice = mean(Practice_T.Correct)*100;
        Edison_Behaviour(i_sub).RT_Practice = mean(Practice_T.Time);
        Edison_Behaviour(i_sub).TotalTime_Practice = TotalTime;
    end
    
    % Pre data
    load(cell2mat(Pre_Sub));Pre_T = struct2table(Pre);
    
    if strcmp(sujet, '47OK') % Nothing in mat, but data in post
        Edison_Behaviour(i_sub).Correct_Pre = nan; Edison_Behaviour(i_sub).RT_Pre = nan;
    else
        Edison_Behaviour(i_sub).Correct_Pre = mean(Pre_T.Correct)*100;
        Edison_Behaviour(i_sub).RT_Pre = mean(Pre_T.Time);
        Edison_Behaviour(i_sub).TotalTime_Pre = TotalTime;
    end
    Edison_Behaviour(i_sub).meanPVT_Pre = meanPVT1;
    
    % For pre solvers
    
    if Edison_Subject.InsightPre(i_sub) ==1
        
        AhaMoment = findchangepts(Pre_T.Time, 'Statistic', 'mean');
        Edison_Behaviour(i_sub).AhaMoment = AhaMoment;
        
        Edison_Behaviour(i_sub).Correct_beforeAha = mean(Pre_T.Correct(1:AhaMoment))*100;
        Edison_Behaviour(i_sub).RT_beforeAha = mean(Pre_T.Time(1:AhaMoment));
        Edison_Behaviour(i_sub).TradeOff_beforeAha = Edison_Behaviour(i_sub).Correct_beforeAha/Edison_Behaviour(i_sub).RT_beforeAha;
        
        Edison_Behaviour(i_sub).Correct_fromAha = mean(Pre_T.Correct(AhaMoment:Nb_Pre_Trial))*100;
        Edison_Behaviour(i_sub).RT_fromAha = mean(Pre_T.Time(AhaMoment:Nb_Pre_Trial));
        Edison_Behaviour(i_sub).TradeOff_fromAha = Edison_Behaviour(i_sub).Correct_fromAha/Edison_Behaviour(i_sub).RT_fromAha;
        
        % Tradeoff for each trial
        for t = 1:length(Pre_T.Correct)
            Edison_Behaviour(i_sub).TradeOff(t) = Pre_T.Correct(t)*100/Pre_T.Time(t);
        end
        
        % Compare tradeoff post start vs post before eureka
        
        if AhaMoment <40 % Not enough trials to compare
            Edison_Behaviour(i_sub).RT_20_Start = nan; Edison_Behaviour(i_sub).RT_20_beforeAha = nan;
            Edison_Behaviour(i_sub).Correct_20_Start = nan; Edison_Behaviour(i_sub).TradeOff_Start =nan;
        else
            Edison_Behaviour(i_sub).RT_20_Start = mean(Pre_T.Time(1:20)); %tradeoff start post
            Edison_Behaviour(i_sub).Correct_20_Start = mean(Pre_T.Correct(1:20))*100;
            Edison_Behaviour(i_sub).TradeOff_Start = Edison_Behaviour(i_sub).Correct_20_Start/Edison_Behaviour(i_sub).RT_20_Start;
            Edison_Behaviour(i_sub).RT_20_beforeAha = mean(Pre_T.Time(AhaMoment-20:AhaMoment)); %tradeoff before aha
        end
        
        % Post analyses on pre solvers
        
        if strcmp(sujet, '25DS') % Pre solver, stopped after pre
            Edison_Behaviour(i_sub).meanPVT_Post = nan;
            Edison_Behaviour(i_sub).Correct_Post = nan;
            Edison_Behaviour(i_sub).RT_Post = nan;
        else
            clear Pre_T;
            load(cell2mat(Post_Sub));Post_T = struct2table(Post);
            Edison_Behaviour(i_sub).meanPVT_Post = meanPVT2;
            Edison_Behaviour(i_sub).Correct_Post = mean(Post_T.Correct(Post_T.Block<10))*100; %before implicit
            Edison_Behaviour(i_sub).RT_Post = mean(Post_T.Time(Post_T.Block<10));
        end
        
        % Post analyses for remaining subjects
        
    elseif strcmp(sujet,'14CD') % Script stopped after insight
        load(cell2mat(Post_Sub));Post_T = struct2table(Post);
        Edison_Behaviour(i_sub).meanPVT_Post = meanPVT2;
        
        nonEmptyAnswer=find(~cellfun(@isempty,Post_T.AnswerString));
        Post_T = Post_T(1:nonEmptyAnswer(end),:);
        Edison_Behaviour(i_sub).Correct_Post = mean(cell2mat(Post_T.Correct(Post_T.Block<10)))*100;
        Edison_Behaviour(i_sub).RT_Post = mean(cell2mat(Post_T.Time(Post_T.Block<10)));
        Edison_Behaviour(i_sub).AhaMoment = 180; AhaMoment = 180; StopTrial = 193;
        Edison_Behaviour(i_sub).Correct_beforeAha = mean(cell2mat(Post_T.Correct(1:AhaMoment-1)))*100;
        Edison_Behaviour(i_sub).RT_beforeAha = mean(cell2mat(Post_T.Time(1:AhaMoment-1)));
        Edison_Behaviour(i_sub).Correct_fromAha = mean(cell2mat(Post_T.Correct(AhaMoment:StopTrial)))*100;
        Edison_Behaviour(i_sub).RT_fromAha = mean(cell2mat(Post_T.Time(AhaMoment:StopTrial)));
        
        Edison_Behaviour(i_sub).TradeOff_fromAha = Edison_Behaviour(i_sub).Correct_fromAha/Edison_Behaviour(i_sub).RT_fromAha;
        Edison_Behaviour(i_sub).TradeOff_beforeAha = Edison_Behaviour(i_sub).Correct_beforeAha/Edison_Behaviour(i_sub).RT_beforeAha;
        
        for t = 1:length(Post_T.Correct)
            Edison_Behaviour(i_sub).TradeOff(t) = cell2mat(Post_T.Correct(t))*100/cell2mat(Post_T.Time(t));
        end
        Edison_Behaviour(i_sub).RT_20_Start = mean(cell2mat(Post_T.Time(1:20)));
        Edison_Behaviour(i_sub).Correct_20_Start = mean(cell2mat(Post_T.Correct(1:20)))*100;
        Edison_Behaviour(i_sub).TradeOff_Start = Edison_Behaviour(i_sub).Correct_20_Start/Edison_Behaviour(i_sub).RT_20_Start;
        Edison_Behaviour(i_sub).RT_20_beforeAha = mean(cell2mat(Post_T.Time(AhaMoment-20:AhaMoment)));
        
    else
        load(cell2mat(Post_Sub));Post_T = struct2table(Post);
        Edison_Behaviour(i_sub).meanPVT_Post = meanPVT2;
        
        if strcmp (sujet,'104MR')
            BeforeImplicit = Post_T.Time(2:Nb_Post_Trial); %first trial: big reaction time, after stable
        elseif strcmp (sujet,'87MD')
            BeforeImplicit = Post_T.Time(1:268); % same with the last two trials
        elseif strcmp (sujet,'71KN')
            BeforeImplicit = Post_T.Time(17:Nb_Post_Trial); %same first trials
        else
            BeforeImplicit = Post_T.Time(Post_T.Block<10);
        end
        
        AhaMoment = findchangepts(BeforeImplicit, 'Statistic', 'mean');
        Edison_Behaviour(i_sub).AhaMoment = AhaMoment;
        Edison_Behaviour(i_sub).Correct_Post = mean(Post_T.Correct(Post_T.Block<10))*100;
        Edison_Behaviour(i_sub).RT_Post = mean(Post_T.Time(Post_T.Block<10));
        
        Edison_Behaviour(i_sub).Correct_fromAha = mean(Post_T.Correct(AhaMoment:Nb_Post_Trial))*100; %270 end of block 9: correct from aha moment to implicite
        Edison_Behaviour(i_sub).RT_fromAha = mean(Post_T.Time(AhaMoment:Nb_Post_Trial));
        Edison_Behaviour(i_sub).TradeOff_fromAha = Edison_Behaviour(i_sub).Correct_fromAha/Edison_Behaviour(i_sub).RT_fromAha;
        
        Edison_Behaviour(i_sub).Correct_beforeAha = mean(Post_T.Correct(1:AhaMoment-1))*100; %270 end of block 9: correct from aha moment to implicite
        Edison_Behaviour(i_sub).RT_beforeAha = mean(Post_T.Time(1:AhaMoment-1));
        Edison_Behaviour(i_sub).TradeOff_beforeAha = Edison_Behaviour(i_sub).Correct_beforeAha/Edison_Behaviour(i_sub).RT_beforeAha;
        
        for t = 1:length(Post_T.Correct)
            Edison_Behaviour(i_sub).TradeOff(t) = Post_T.Correct(t)*100/Post_T.Time(t);
        end
        if AhaMoment <40
            Edison_Behaviour(i_sub).RT_20_Start = nan;
            Edison_Behaviour(i_sub).RT_20_beforeAha = nan;
            Edison_Behaviour(i_sub).TradeOff_Start= nan;
            Edison_Behaviour(i_sub).Correct_20_Start = nan;
        else
            Edison_Behaviour(i_sub).RT_20_Start = mean(Post_T.Time(1:20));
            Edison_Behaviour(i_sub).Correct_20_Start = mean(Post_T.Correct(1:20))*100;
            Edison_Behaviour(i_sub).TradeOff_Start = Edison_Behaviour(i_sub).Correct_20_Start/Edison_Behaviour(i_sub).RT_20_Start;
            Edison_Behaviour(i_sub).RT_20_beforeAha = mean(Post_T.Time(AhaMoment-20:AhaMoment));
        end
    end
    
    %------------------------------------------
    % Implicit analyses
    %--------------------------------------------
    
    if strcmp(sujet,'14CD') || strcmp(sujet, '25DS') % Script stopped after insight or pre solver (without post data)
        Edison_Behaviour(i_sub).Correct_Implicit = NaN; Edison_Behaviour(i_sub).RT_Post_Implicit = NaN;
        Edison_Behaviour(i_sub).Correct_Implicit_Rule = nan; Edison_Behaviour(i_sub).RT_Implicit_Rule =nan;
        Edison_Behaviour(i_sub).Correct_Implicit_BreakRule =nan; Edison_Behaviour(i_sub).RT_Implicit_BreakRule = nan;
        Edison_Behaviour(i_sub).TradeOff_Implicit_Rule =nan; Edison_Behaviour(i_sub).TradeOff_Implicit_BreakRule =nan;
        Edison_Behaviour(i_sub).TotalTime_Post = NaN;
        
        
    else
        
        Edison_Behaviour(i_sub).Correct_Implicit = mean(Post_T.Correct(Post_T.Block>9))*100;
        Edison_Behaviour(i_sub).RT_Post_Implicit = mean(Post_T.Time(Post_T.Block>9));
        
        Edison_Behaviour(i_sub).Correct_Implicit_Rule = mean(Post_T.Correct(Post_T.Block>9 & strcmp(Post_T.Type,'Rule')))*100;
        Edison_Behaviour(i_sub).RT_Implicit_Rule = mean(Post_T.Time(Post_T.Block>9 & strcmp(Post_T.Type,'Rule')));
        Edison_Behaviour(i_sub).TradeOff_Implicit_Rule =  Edison_Behaviour(i_sub).Correct_Implicit_Rule/Edison_Behaviour(i_sub).RT_Implicit_Rule;
        
        Edison_Behaviour(i_sub).Correct_Implicit_BreakRule = mean(Post_T.Correct(Post_T.Block>9 & strcmp(Post_T.Type,'BreakRule')))*100;
        Edison_Behaviour(i_sub).RT_Implicit_BreakRule = mean(Post_T.Time(Post_T.Block>9 & strcmp(Post_T.Type,'BreakRule')));
        Edison_Behaviour(i_sub).TradeOff_Implicit_BreakRule =  Edison_Behaviour(i_sub).Correct_Implicit_BreakRule/Edison_Behaviour(i_sub).RT_Implicit_BreakRule;
        
        Edison_Behaviour(i_sub).TotalTime_Post = TotalTime;
       
    end
    
end

% Save
Edison_Behaviour = struct2table(Edison_Behaviour);
Data = [Edison_Subject, Edison_EEG, Edison_Behaviour];

save(strcat(SavePath,'Edison_Behaviour'), 'Edison_Behaviour');
save(strcat(SavePath,'Edison_Subject'), 'Edison_Subject');
save(strcat(SavePath,'All_Data'), 'Data');



