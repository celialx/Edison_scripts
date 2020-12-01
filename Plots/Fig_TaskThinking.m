Init;

figure

decalage = 0;

N = 4;
X = linspace(0,pi*3,1000);
C = linspecer(N);

stackedTask = [Nb_TaskThinking-Nb_TaskThinking_Solver,Nb_TaskThinking_Solver];
stackedHypna = [Nb_Hypna-Nb_Hypna_Solver, Nb_Hypna_Solver];
stackedHypnaTask = [Nb_HypnaTask-Nb_HynaTask_Solver, Nb_HynaTask_Solver];

figure
bh = bar([1 2 3], [stackedTask; stackedHypna; stackedHypnaTask], 'stacked', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2);

for i =1:length(bh)
    set(bh(i), 'FaceColor', ColorSolver2(i,:))
end

leg = legend({'Non Solvers', 'Solvers'});
legend boxoff

Labels = {'TaskThinking', 'Hypna', 'HypnaTask'};
TextLabel = 'Repartition of subjects (%)';

%
[h,p] = prop_test([Nb_TaskThinking_Solver, Nb_TaskThinking-Nb_TaskThinking_Solver], [Nb_Solvers_Post,Nb_NonSolvers_Post], true);
higherBar = 0; sigstar([1 1.1], p, [], higherBar)

[h,p] = prop_test([Nb_Hypna_Solver,Nb_Hypna-Nb_Hypna_Solver], [Nb_Solvers_Post,Nb_NonSolvers_Post], true);
higherBar = 0; sigstar([2 2.1], p, [], higherBar)

[h,p] = prop_test([Nb_HynaTask_Solver,Nb_HypnaTask-Nb_HynaTask_Solver], [Nb_Solvers_Post,Nb_NonSolvers_Post], true);
higherBar = 0; sigstar([3 3.1], p, [], higherBar)

Design;


allChildren = get(gca, 'Children');
displayNames = get(allChildren, 'DisplayName');
delete(allChildren(strcmp(displayNames, 'data1'))); delete(allChildren(strcmp(displayNames, 'data2')))
delete(allChildren(strcmp(displayNames, 'data3')))

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_TaskThinking.jpg')],'-djpeg','-r600');
end
%% Hypna vs insight
Init;
decalage = 0; higherBar = 0;

HypnaTask_WakeSolvers = nnz(HypnaTask_Solver & WakeAll)/Nb_Wake_Solver*100;
HypnaTask_N1Solvers =  nnz(HypnaTask_Solver & N1All)/Nb_N1_Solver*100;

Hypna_WakeSolvers = nnz(Hypna_Solver & WakeAll)/Nb_Wake_Solver*100;
Hypna_N1Solvers = nnz(Hypna_Solver & N1All)/Nb_N1_Solver*100;

Task_Thinking_WakeSolvers = nnz(TaskThinking_Solver & WakeAll)/Nb_Wake_Solver*100;
TaskThinking_N1Solvers = nnz(TaskThinking_Solver & N1All)/Nb_N1_Solver*100;

Var = [HypnaTask_WakeSolvers,HypnaTask_N1Solvers,Hypna_WakeSolvers, Hypna_N1Solvers, Task_Thinking_WakeSolvers, TaskThinking_N1Solvers];

% for i = 1:6
%     
%     bar(i, Var(i));
%     hold on
% end

% TextLabel = 'Percentage of solvers';
% Labels = 'Trials'; 
% Design;

