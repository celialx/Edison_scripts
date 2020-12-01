
Init;
higherBar= 0; decalage = 1;

Var = {Data.TradeOff_Implicit_Rule(AllSol), Data.TradeOff_Implicit_BreakRule(AllSol),...
    Data.TradeOff_Implicit_Rule(AllNonSol), Data.TradeOff_Implicit_BreakRule(AllNonSol)}

%     Var = {Impli_Sol_Correct_Rule, Impli_Sol_Correct_BreakRule, Impli_NonSol_Correct_Rule, Impli_NonSol_Correct_BreakRule, Impli_Sol_RT_Rule...
%         , Impli_Sol_RT_BreakRule, Impli_NonSol_RT_Rule, Impli_NonSol_RT_BreakRule};

figure

for i =1:2
    
    subplot(1,2,i);
    N = length(Var{i+1}); X = linspace(0,pi*3,1000); C = linspecer(N);
    
    if i ==1
        title('Solvers')
        
        for ii=1:N
            plot([Var{i}(ii), Var{i+1}(ii)]','LineWidth', 2, 'color', C(ii,:), 'linestyle', '--')
            hold on;
        end
        
        [h,p] = ttest(Var{i}, Var{i+1});
        sigstar([1 2], p, [],higherBar);
        
    elseif i ==2
        title('Non Solvers')
        
        for ii=1:N
            plot([Var{i+1}(ii), Var{i+2}(ii)]','LineWidth', 2, 'color', C(ii,:), 'linestyle', '--')
            hold on;
        end
        
        [h,p] = ttest(Var{i+1}, Var{i+2});
        sigstar([1 2], p, [],higherBar);
    end
    
    TextLabel = 'Speed-Accuracy Trade-off';
    Labels = {'Rule', '\0Rule'};
    set(gca, 'xlim', [0.5 2.5]);
    
    Design;
    
    
end

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Implicit.jpg')],'-djpeg','-r600');
end

%% Figure plotting the difference between both trials of implicit
%
Init;
decalage =1; higherBar =0;

figure
diff_Sol_Correct(isnan(diff_Sol_Correct)) = [];
diff_Sol_RT(isnan(diff_Sol_RT)) = [];

Var = {diff_Sol_Correct, diff_NonSol_Correct, diff_Sol_RT, diff_NonSol_RT};
cte = 0;
for i =1:2
    if i ==1
        cte= cte+1;
    else
        cte= cte+2;
    end
    
    subplot(1,2,i)
    x = [ones(length(Var{cte}),1); 1+ones(length(Var{cte+1}),1)];
    y = [Var{cte}; Var{cte+1}];
    h = beeswarm(x, y, 'sort_style', 'up', 'overlay_style', 'ci', 'dot_size', 1.3, 'colormap', ColorSolver, 'MarkerFacealpha', 0.5,'MarkerEdgeColor','k');
    hold off
    Labels = LabelsSol;
    if i==1
        TextLabel = ' Accuracy (%) : Rule - \0Rule';
    else
        TextLabel = 'Solving Time (secs) : Rule - \0Rule';
    end
    
    yline(0,'--', 'LineWidth', 2.5,'Color',Basiccolor);
    [h,p] = ttest2(Var{cte}, Var{cte+1});
    sigstar([1 2], p ,[], higherBar);
    
    Design;
    
end


if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_diff_Impli_correct.jpg')],'-djpeg','-r600');
end

%% Figure plotting the difference between both trials of implicit


Init;
decalage =1; higherBar =0;

figure

Var = {diff_NonSol_Correct_Wake; diff_NonSol_Correct_N1; diff_NonSol_Correct_N2;...
    diff_NonSol_RT_Wake; diff_NonSol_RT_N1; diff_NonSol_RT_N2};
cte = 0;
for i =1:2
    if i ==1
        cte= cte+1;
    else
        cte= cte+3;
    end
    
    subplot(1,2,i)
    x = [ones(length(Var{cte}),1); 1+ones(length(Var{cte+1}),1); 2+ones(length(Var{cte+2}),1)];
    y = [Var{cte}; Var{cte+1}; Var{cte+2}];
    h = beeswarm(x, y, 'sort_style', 'up', 'overlay_style', 'ci', 'dot_size', 1.3, 'colormap', Color3N2, 'MarkerFacealpha', 0.5,'MarkerEdgeColor','k');
    hold off
    Labels = Labels3N2;
    if i==1
        TextLabel = ' Accuracy (%) : Rule - \0Rule';
    else
        TextLabel = 'Solving Time (secs) : Rule - \0Rule';
    end
    
    yline(0,'--', 'LineWidth', 2.5,'Color',Basiccolor);
    p = anovan(y,x, 'display', 'off');
    sigstar([1 3], p ,[], higherBar);
    
    Design;
end


if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_diff_Impli_Correct.jpg')],'-djpeg','-r600');
end


%% Same for tradeoff
%
Init;
decalage =0; higherBar =0;

figure

diff_TradeOff_Sol = minus(Data.TradeOff_Implicit_Rule(AllSol),  Data.TradeOff_Implicit_BreakRule(AllSol));
diff_TradeOff_NonSol = minus(Data.TradeOff_Implicit_Rule(AllNonSol),  Data.TradeOff_Implicit_BreakRule(AllNonSol))

diff_TradeOff_Sol(isnan(diff_TradeOff_Sol)) = [];
diff_TradeOff_NonSol(isnan(diff_TradeOff_NonSol)) = [];

Var = {diff_TradeOff_Sol, diff_TradeOff_NonSol};

x = [ones(length(Var{1}),1); 1+ones(length(Var{2}),1)];
y = [Var{1}; Var{2}];
h = beeswarm(x, y, 'sort_style', 'up', 'overlay_style', 'ci', 'dot_size',1.3, 'colormap', ColorSolver, 'MarkerFacealpha', 0.5,'MarkerEdgeColor','k');
Labels = LabelsSol;
TextLabel = 'Speed-Accuracy Trade-off';

yline(0,'--', 'LineWidth', 2.5,'Color',Basiccolor);
[h,p] = ttest2(Var{1}, Var{2});
sigstar([1 2], p ,[], higherBar);


Design;

figure
diff_TradeOff_NonSol_N1 = minus(Data.TradeOff_Implicit_Rule(AllNonSol & N1All),  Data.TradeOff_Implicit_BreakRule(AllNonSol & N1All)); diff_TradeOff_NonSol_N1(isnan(diff_TradeOff_NonSol_N1)) = [];
diff_TradeOff_NonSol_Wake = minus(Data.TradeOff_Implicit_Rule(AllNonSol & WakeAll),  Data.TradeOff_Implicit_BreakRule(AllNonSol & WakeAll)); diff_TradeOff_NonSol_Wake(isnan(diff_TradeOff_NonSol_Wake)) = [];
diff_TradeOff_NonSol_N2 = minus(Data.TradeOff_Implicit_Rule(AllNonSol & N2All),  Data.TradeOff_Implicit_BreakRule(AllNonSol & N2All)); diff_TradeOff_NonSol_N2(isnan(diff_TradeOff_NonSol_N2)) = [];

Var = {diff_TradeOff_NonSol_Wake; diff_TradeOff_NonSol_N1; diff_TradeOff_NonSol_N2};

x = [ones(length(Var{1}),1); 1+ones(length(Var{2}),1); 2+ones(length(Var{3}),1)];
y = [Var{1}; Var{2}; Var{3}];
h = beeswarm(x, y, 'sort_style', 'up', 'overlay_style', 'ci', 'dot_size', 1.3, 'colormap', Color3N2, 'MarkerFacealpha', 0.5,'MarkerEdgeColor','k');
hold off
Labels = Labels3N2;
yline(0,'--', 'LineWidth', 2.5,'Color',Basiccolor);
p = anovan(y,x, 'display','off');
sigstar([1 3], p ,[], higherBar);
Design;


if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_diff_Impli_TradeOff.jpg')],'-djpeg','-r600');
end

