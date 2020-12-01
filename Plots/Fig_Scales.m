% FIg scales

% For numerical data

Init;
fig = figure
Var = [Data.Motivation_PostStart, Data.Focus_PostStart, Data.Boredom_PostStart, Data.Sleepiness_PostStart];

Factor_label = {'Motivation'; 'Concentration'; 'Boredom'; 'Sleepiness'};

decalage = 0;

for i =1:size(Var,2)
    
    Factor = Var(:,i);
    Var1 = Factor(Wake); Var2 = Factor(N1);
    
    subplot(2,2,i);
    
    violinPlot(Var1,'histOri', 'left', 'widthDiv', [2 1], 'showMM', 2, ...
        'color', Color2(1,:));
    
    hold on;
    
    violinPlot(Var2, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 2, ...
        'color', Color2(2,:));
    
    TextLabel = cell2mat(Factor_label(i));
    Labels = 'Trials';
    set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'Wake', 'N1'}, 'xlim', [0.2 1.8]);
    Design;
    
    [~, pval] = ttest2(Var1, Var2);
    mysigstar(gca, xticks, 3/4*ylim(2), pval);
    
end

Design;
% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off');

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Scales.jpg')],'-djpeg','-r600');
end
