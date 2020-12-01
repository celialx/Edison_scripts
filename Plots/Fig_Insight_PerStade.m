%% Figure 1 : Percentage of Insight per Group

Init;

figure;
subplot(1,2,1)
decalage = 1; higherBar = 0;

Var = [Pre_Solver_Percent, Post_Solver_Wake_Percent, Post_Solver_N1_Percent, Post_Solver_N2_Percent];

for i = 1:4
    bar(i, Var(i), 'FaceColor', Color4(i,:), 'EdgeColor', Edge4(i,:), 'LineWidth', 2);
    hold on
end

TextLabel = 'Insight (%)';
Labels = Labels4;

Design;

y1 = set(gca, 'ylim',[0 100]);
yticks([0:20:100]);

% Stat: contingency table (fischer test)
% v =1: comparison of pre with remaining groups
% v =2: comparison of wake with remaining groups

for v = 1:2
    
    if v ==1
        Var_cte_solver = Nb_Pre_Solver;
        Var_cte = Nb_Sub;
        
        Var_diff_solver = [Nb_Wake_Solver, Nb_N1_Solver, Nb_N2_Solver];
        Var_diff = [Nb_Wake, Nb_N1, Nb_N2];
        ColNames = {'Solvers', 'NonSolvers'};
        Name_cte = 'Pre';
        Name = {'Wake', 'N1', 'N2'};
        n =1:3;
        
    elseif v ==2
        Var_cte_solver = Nb_Wake_Solver;
        Var_cte = Nb_Wake;
        
        Var_diff_solver = [Nb_N1_Solver, Nb_N2_Solver];
        Var_diff = [Nb_N1, Nb_N2];
        ColNames = {'Solvers', 'NonSolvers'};
        Name_cte = 'Wake';
        Name = {'N1', 'N2'};
        n =1:2;
    end
    
    for i =1:length(n)
        
        table_cont = table([Var_cte_solver; Var_diff_solver(i)], [Var_cte - Var_cte_solver; Var_diff(i) - Var_diff_solver(i)],...
            'VariableNames', ColNames, 'RowNames', {Name_cte, Name{i}})
        
        [h,p] = fishertest(table_cont);
        if v ==1
            sigstar([1,i+1], p, [], higherBar);
        elseif v==2
            sigstar([2,i+2], p, [], higherBar);
        end
    end
end

% N1 vs N2

table_N1N2 = table([Nb_N1_Solver; Nb_N2_Solver],[Nb_N1 - Nb_N1_Solver ; Nb_N2 - Nb_N2_Solver],...
    'VariableNames', ColNames, 'RowNames', {'Post N1', 'Post N2'})

[h,p] = fishertest(table_N1N2);
sigstar([3,4], p, [], higherBar);

% Same with stacked bar (distribution of solvers vs nonsol in each group)


subplot(1,2,2)
stacked_Wake = [Nb_WakeAll - (Nb_Pre_Wake_Solver + Nb_Wake_Solver), Nb_Pre_Wake_Solver, Nb_Wake_Solver];
stacked_N1 = [Nb_N1All - (Nb_Pre_N1_Solver + Nb_N1_Solver), Nb_Pre_N1_Solver, Nb_N1_Solver];
stacked_N2 = [Nb_N2All - (Nb_Pre_N2_Solver + Nb_N2_Solver), Nb_Pre_N2_Solver, Nb_N2_Solver];

c = cbrewer('seq', 'OrRd',3);

bh = bar([1 2 3],[stacked_Wake; stacked_N1; stacked_N2], 'stacked', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 2.5);

for i =1:length(bh)
    set(bh(i), 'FaceColor', c(i,:));
end

TextLabel = 'Repartition of subjects'; Labels = Labels3N2;
legend(' Non Solvers', ' Solvers{\it Pre}', ' Solvers{\it Post}'); legend boxoff;
Design;

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Insight.jpg')],'-djpeg','-r600');
end