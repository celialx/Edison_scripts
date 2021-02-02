clear all; close all; clc;

%-----------------------------------------------------------------
% Set the paths where the data are stored
%------------------------------------------------------------------

FigPath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Figures\';
TablePath = 'C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Edison_Tables\';
load([TablePath, 'All_Data']);
path_local ='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\LocalSleep\';
path_LSCPtools='D:\MATLAB\Toolbox\LSCPtools\';
path_localsleep='D:\MATLAB\Toolbox\wanderIM\localsleep\';
addpath(path_localsleep);
addpath(genpath(path_LSCPtools));

%If you desire to print figure, set this to 1
SaveFig = 0;

%-----------------------------------------------------------------
% Initialisation
%------------------------------------------------------------------

Basiccolor =  [0.2 0.2 0.2];
orange = [1 0.6 0.6]; orangeedge = [1 0.6 0.4];
verteau = [0.2 1 0.8]; verteauedge = [0.4 0.8 0.6]; % Solver color
verteau2 = [0 0.8 0.4];verteauedge2 = [0 0.6 0.4];
red = [1 0 0.2]; rededge = [1 0 0]; % NonSol color

precolor = [1 0.69 0.4]; preedge = [1 0.41 0.16];
wakecolor = [0 1 0.45]; wakeedge = [0.24 0.70 0.15];
N1color =  [0 1 1]; N1edge = [0.11 0.50 0.50];
N2color = [0 0.5 0.5]; N2edge = [0.16 0.23 0.27];

Color2 = [wakecolor; N1color];Edge2 = [wakeedge; N1edge];
Color3 = [precolor; wakecolor; N1color]; Edge3 = [preedge; wakeedge; N1edge];
Color3N2 = [wakecolor; N1color; N2color]; Edge3N2 = [wakeedge; N1edge; N2edge];
Color4 = [precolor; wakecolor; N1color; N2color]; Edge4 = [preedge; wakeedge; N1edge; N2edge];
ColorImplicit = [verteau; wakecolor; N2color; red];

ColorN1 = [N1color; N2color];
ColorSolver = [verteau; red]; ColorNonSol = [red; verteau; N1color ];
ColorSolver2 = [red; verteau];

% Labels for x axes
LabelsBall = {'AASM'; 'BERN'}; LabelsSol = {'Solvers'; 'Non Solvers'};
Labels2 = {'Wake'; 'N1'}; Labels3 = {'\it{Pre}'; 'Wake'; 'N1'}; 
Labels3N2 = {'Wake'; 'N1'; 'N2'}; Labels4 = {'\it{Pre}'; 'Wake'; 'N1'; 'N2'};

%-----------------------------------------------------------------
% Groups
%------------------------------------------------------------------

Nb_Sub = length(Data.Code); % All subjects

% Groups (all subjects with pre)
WakeAll = Data.HypnoEdisonN1==0 & Data.HypnoEdisonN2 ==0; Nb_WakeAll = nnz(WakeAll);
N1All = Data.HypnoEdisonN1>0 & Data.HypnoEdisonN2 ==0; Nb_N1All = nnz(N1All);
N2All  = Data.HypnoEdisonN2>0; Nb_N2All = nnz(N2All);
N1_N2All =  N1All | N2All; Nb_N1_N2All = nnz(N1_N2All);

% Pre solvers
Pre_Solver = Data.InsightPre==1; Nb_Pre_Solver = nnz(Pre_Solver); 
Pre_Wake_Solver = Pre_Solver & WakeAll; Nb_Pre_Wake_Solver = nnz(Pre_Wake_Solver);
Pre_N1_Solver = Pre_Solver & N1All; Nb_Pre_N1_Solver = nnz(Pre_N1_Solver);
Pre_N2_Solver = Pre_Solver & N2All; Nb_Pre_N2_Solver = nnz(Pre_N2_Solver);

Pre_Solver_Percent = (Nb_Pre_Solver/Nb_Sub)*100; 
Pre_Solver_Wake_Percent = (Nb_Pre_Wake_Solver/Nb_WakeAll)*100;
Pre_Solver_N1_Percent = (Nb_Pre_N1_Solver/Nb_N1All)*100;
Pre_Solver_N2_Percent = (Nb_Pre_N2_Solver/Nb_N2All)*100;

% Post subjects without pre subjects
Wake = WakeAll & Data.InsightPre==0; Nb_Wake = nnz(Wake);
N1 = N1All & Data.InsightPre==0; Nb_N1 =nnz(N1); % Pure N1, without N2 intrusions
N2 = N2All & Data.InsightPre==0; Nb_N2 = nnz(N2); % All N2 subjects (N1+N2 & N2 only)

% Post solvers (without taking into account pre subjects)
Wake_Solver = Data.InsightPost(Wake); Nb_Wake_Solver = nnz(Wake_Solver); Post_Solver_Wake_Percent = (Nb_Wake_Solver/Nb_Wake)*100;
N1_Solver = Data.InsightPost(N1); Nb_N1_Solver = nnz(N1_Solver); Post_Solver_N1_Percent = (Nb_N1_Solver/Nb_N1)*100;
N2_Solver =Data.InsightPost(N2); Nb_N2_Solver = nnz(N2_Solver); Post_Solver_N2_Percent = (Nb_N2_Solver/Nb_N2)*100;

% Data on sol vs non sol
AllSol = Data.InsightPre ==1 | Data.InsightPost==1;
SolPost = Data.InsightPre ==0 & Data.InsightPost==1; Nb_Solvers_Post = nnz(SolPost);
AllNonSol = Data.InsightPre ==0 & Data.InsightPost==0; Nb_NonSolvers_Post = nnz(AllNonSol);
AllSub = AllNonSol | AllSol;
AllPost = SolPost | AllNonSol;
% Same without N2 subjects
AllSolN2 = AllSol & Data.HypnoEdisonN2 ==0;
SolPostN2 = SolPost & Data.HypnoEdisonN2 ==0;
AllNonSolN2 = AllNonSol & Data.HypnoEdisonN2 ==0;

%-----------------------------------------------------------------
% Analyses on bottle
%------------------------------------------------------------------

Nb_Ball = sum(Data.Fall); % Nb Total of Ball fell (with the repetitions per subject)
Nb_Ball_Single = nnz(Data.Fall>0); % Nb Total of Ball fell (without repetitions, limited to 1 per subject)

% At which stade the bottle dropped (looked 2 secs before bottle dropped)
stadeBall_Multiple = cell2mat(Data.stadeBall(Data.Fall>1)'); % ball dropped several times
stadeBall_One = cell2mat(Data.stadeBall(Data.Fall==1)'); % ball dropped one time
stadeBall_All = [stadeBall_Multiple, stadeBall_One]'; %altogether

N1BfBall = nnz(stadeBall_All==1)/Nb_Ball*100;
WakeBfBall = nnz(stadeBall_All==0)/Nb_Ball*100;

% After how long did it fall?
timeN1bfBall = horzcat(Data.timeN1bfBall{:}); alltimeN1bfBall = timeN1bfBall;
timeN1bfBall(isnan(timeN1bfBall)) =[];
%-----------------------------------------------------------------
% Hypnagogia Bottle vs No bottle
%------------------------------------------------------------------

NoBottle = Data.Bottle==0; Nb_NoBottle = nnz(NoBottle);
HypnaGeneral = nnz(Data.Hypnagogia ==1 & NoBottle); 

Bottle = Data.Bottle ==1; Nb_Bottle = nnz(Bottle);
HypnaBottle = nnz(Data.HypnaBottle ==1 & Bottle); 

%-----------------------------------------------------------------
% Task thinking, hypnagogia on task
%------------------------------------------------------------------

% without pre solvers
HypnaTask = Data.HypnaTask ==1 & Data.InsightPre ==0;
Nb_HypnaTask = nnz(HypnaTask); Percent_HypnaTask =Nb_HypnaTask/nnz(AllPost)*100;

HypnaTask_Solver = HypnaTask & Data.InsightPost ==1; Nb_HynaTask_Solver = nnz(HypnaTask_Solver);
Nb_HypnaTask_Solver_Wake = nnz(HypnaTask_Solver & WakeAll);

Hypna = Data.Hypnagogia ==1 & Data.InsightPre ==0; 
Nb_Hypna = nnz(Hypna); Percent_Hypna = Nb_Hypna/ nnz(AllPost)*100;

Hypna_Solver = Hypna & Data.InsightPost ==1;Nb_Hypna_Solver = nnz(Hypna_Solver);
Hypna_Solver_Wake = nnz(Hypna_Solver & WakeAll);

TaskThinking = Data.TaskThinking ==1 & Data.InsightPre ==0;
Nb_TaskThinking = nnz(TaskThinking); Percent_TaskThinking =Nb_TaskThinking / nnz(AllPost)*100;
TaskThinking_Solver = TaskThinking & Data.InsightPost ==1; Nb_TaskThinking_Solver = nnz(TaskThinking_Solver);

%-----------------------------------------------------------------
% Aha Moment
%------------------------------------------------------------------

% At which trial does the Eureka occur?
AhaMoment = Data.AhaMoment(SolPost);% (restreined here to post solvers)

% Is there a difference btw vigilance states?
AhaMoment_Wake = Data.AhaMoment(Data.InsightPost==1 & WakeAll);
AhaMoment_N1 = Data.AhaMoment(Data.InsightPost==1 & N1All);
AhaMoment_N2 = Data.AhaMoment(Data.InsightPost==1 & N2All);

%-----------------------------------------------------------------
% Analyses on implicit task
%------------------------------------------------------------------

%----------------------------------------
% For Solvers subjects (pre or post)
%----------------------------------------

% Perf on rule trials
Impli_Sol_Correct_Rule = Data.Correct_Implicit_Rule(AllSol); 
Impli_Sol_RT_Rule = Data.RT_Implicit_Rule(AllSol); 

% Perf on break trials
Impli_Sol_Correct_BreakRule = Data.Correct_Implicit_BreakRule(AllSol); 
Impli_Sol_RT_BreakRule = Data.RT_Implicit_BreakRule(AllSol); 

% Difference between correct and break trials
diff_Sol_RT = minus(Impli_Sol_RT_Rule, Impli_Sol_RT_BreakRule);
diff_Sol_Correct = minus(Impli_Sol_Correct_Rule, Impli_Sol_Correct_BreakRule);

%----------------------------------------
% For NonSolvers subjects
%----------------------------------------

% Perf rule trials
Impli_NonSol_Correct_Rule = Data.Correct_Implicit_Rule(AllNonSol);
Impli_NonSol_RT_Rule = Data.RT_Implicit_Rule(AllNonSol);
% Subgroups
Impli_NonSol_Correct_Rule_Wake = Data.Correct_Implicit_Rule(AllNonSol & WakeAll);
Impli_NonSol_RT_Rule_Wake = Data.RT_Implicit_Rule(AllNonSol & WakeAll);

Impli_NonSol_Correct_Rule_N1 = Data.Correct_Implicit_Rule(AllNonSol & N1All);
Impli_NonSol_RT_Rule_N1 = Data.RT_Implicit_Rule(AllNonSol & N1All);

Impli_NonSol_Correct_Rule_N2 = Data.Correct_Implicit_Rule(AllNonSol & N2All);
Impli_NonSol_RT_Rule_N2 = Data.RT_Implicit_Rule(AllNonSol & N2All);

% Perf on break trials
Impli_NonSol_Correct_BreakRule = Data.Correct_Implicit_BreakRule(AllNonSol);
Impli_NonSol_RT_BreakRule = Data.RT_Implicit_BreakRule(AllNonSol);
%Subgroups
Impli_NonSol_Correct_BreakRule_Wake = Data.Correct_Implicit_BreakRule(AllNonSol & WakeAll);
Impli_NonSol_RT_BreakRule_Wake = Data.RT_Implicit_BreakRule(AllNonSol & WakeAll);

Impli_NonSol_Correct_BreakRule_N1 = Data.Correct_Implicit_BreakRule(AllNonSol & N1All);
Impli_NonSol_RT_BreakRule_N1 = Data.RT_Implicit_BreakRule(AllNonSol & N1All);

Impli_NonSol_Correct_BreakRule_N2 = Data.Correct_Implicit_BreakRule(AllNonSol & N2All);
Impli_NonSol_RT_BreakRule_N2 = Data.RT_Implicit_BreakRule(AllNonSol & N2All);

% Difference between correct and break trials
diff_NonSol_Correct = minus(Impli_NonSol_Correct_Rule, Impli_NonSol_Correct_BreakRule);
diff_NonSol_RT= minus(Impli_NonSol_RT_Rule, Impli_NonSol_RT_BreakRule);

diff_NonSol_Correct_Wake = minus(Impli_NonSol_Correct_Rule_Wake, Impli_NonSol_Correct_BreakRule_Wake);
diff_NonSol_RT_Wake= minus(Impli_NonSol_RT_Rule_Wake, Impli_NonSol_RT_BreakRule_Wake);

diff_NonSol_Correct_N1= minus(Impli_NonSol_Correct_Rule_N1, Impli_NonSol_Correct_BreakRule_N1);
diff_NonSol_RT_N1= minus(Impli_NonSol_RT_Rule_N1, Impli_NonSol_RT_BreakRule_N1);

diff_NonSol_Correct_N2= minus(Impli_NonSol_Correct_Rule_N2, Impli_NonSol_Correct_BreakRule_N2);
diff_NonSol_RT_N2= minus(Impli_NonSol_RT_Rule_N2, Impli_NonSol_RT_BreakRule_N2);
