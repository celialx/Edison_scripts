% FIg scales

% For numerical data

Init;
fig = figure
Var = [Data.Motivation_PostStart, Data.Focus_PostStart, Data.Boredom_PostStart, Data.Sleepiness_PostStart];

Factor_label = {'Motivation'; 'Concentration'; 'Boredom'; 'Sleepiness'};

decalage = 0; higherBar = 0;

for i =1:size(Var,2)
    
    Factor = Var(:,i);
    Var1 = Factor(Wake); Var2 = Factor(N1); Var3 = Factor(N2);
    
    subplot(2,2,i);
    
    bplot(Var1, 1, 'color',Color3N2(1,:));
%     violinPlot(Var1,'histOri', 'left', 'widthDiv', [2.5 1], 'showMM', 2, ...
%         'color', Color3N2(1,:));
    
    hold on;
        bplot(Var2, 2, 'color',Color3N2(2,:));

%     violinPlot(Var2, 'histOri', 'right', 'widthDiv', [2.5 2], 'showMM', 2, ...
%         'color', Color3N2(2,:));
    
    hold on 
            bplot(Var3, 3, 'color',Color3N2(3,:));

%     violinPlot(Var3, 'histOri', 'right', 'widthDiv', [2.5 2.5], 'showMM', 2, ...
%         'color', Color3N2(3,:));
    
    TextLabel = cell2mat(Factor_label(i));
    Labels = Labels3N2;
%     set(gca, 'xtick', [0.6 1.4 2.2], 'xticklabel', {'Wake', 'N1', 'N2'}, 'xlim', [0.2 1.8]);
    Design;
    
    data = [Var1; Var2; Var3];
    groups = [zeros(1,length(Var1))'; ones(1, length(Var2))'; 1+ones(1,length(Var3))'];
    p = kruskalwallis(data, groups, 'off');
    sigstar([1 3],p, [], higherBar)
%     mysigstar(gca, xticks, 3/4*ylim(2), pval);
    
end

Design;
% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off');

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Scales.jpg')],'-djpeg','-r600');
end

% For kruskal
% data = [Data.Boredom_PostStart(Wake);Data.Boredom_PostStart(N1); Data.Boredom_PostStart(N2)]
% groups = [zeros(1,length(Data.Boredom_PostStart(Wake)))'; ones(1,length(Data.Boredom_PostStart(N1)))'; 1+ones(1,length(Data.Boredom_PostStart(N2)))'];
% [p, ANOVATAB, STATS] = kruskalwallis(data, groups)