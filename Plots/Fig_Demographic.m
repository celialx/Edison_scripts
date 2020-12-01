%-----------------------------------------------------------------
% Figure - Descriptive statistics
%------------------------------------------------------------------
Init;

% For numerical data
Factor_num = [Data.Age, Data.NC, Data.Epworth];
Factor_label = {'Age'; 'Educational level'; 'Epworth Score'};

figure
decalage =1;

for i = 1:size(Factor_num,2)
    
    fig =  subplot(1,3,i)
    
    Factor = Factor_num(:,i);
    Var1 = Factor(WakeAll); Var2= Factor(N1All);
    violinPlot(Var1, 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 2, ...
        'color', Color2(1,:));
    
    hold on;
    
    violinPlot(Var2, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 2, ...
        'color', Color2(2,:));
    
    Labels = 'Trials';
    TextLabel = cell2mat(Factor_label(i));
    set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'Wake', 'N1'}, 'xlim', [0.2 1.8]);
    Design;
    
    
    % Quantitative data
    [~, pval] = ttest2(Var1, Var2);
    mysigstar(gca, xticks, 3/4*ylim(2), pval);
    
end

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Demographic1.jpg')],'-djpeg','-r600');
end

%% For categorical data
% Stat pour données ordinales (un peu, beaucoup... : kruskal_wallis ordinal
% Stat pour données qualitatives :chi square (prop_test)
Init;

Factor_cate = [Data.Gender, Data.Laterality, Data.FieldClass, Data.DreamRecall];
Factor_label2 = {'Gender'; 'Laterality'; 'Field of Expertise'; 'Dream Recall Frequency'};

fig = figure
decalage = 1;

for i = 1:size(Factor_cate,2)
    
    a = numSubplots(size(Factor_cate,2));
    subplot(a(1),a(2),i);
    
    map = Color3N2; colormap(map);
    
    % Define variables
    Factor = Factor_cate(:,i);
    you.Wake = Factor(WakeAll); you.N1 = Factor(N1All); you.N2 = Factor(N2All);
    
    [a,b,c] = nhist(you, 'p','Proportion', 'fsize', 18, 'linewidth', 4, 'smooth', 'color','colormap')
    hold off;
    
    groups = [zeros(1,length(b.Wake)-2) ones(1,length(b.N1)-2) 1+ones(1,length(b.N2)-2)];
    data = [b.Wake(1:end-2) b.N1(1:end-2) b.N2(1:end-2)];
    
    if i == 3 || i == 4 || i ==5 % Ordinals more than two groups
        pkrus(i) = kruskalwallis(data, groups, 'off');
        %
        %     elseif i ==1 || i ==2 % Qualitatives nominales: two groups
        %         [h(i) pchi(i)] = prop_test([b.Wake(1) b.N1(1) b.N2(1)], [b.Wake(1) + b.Wake(2), b.N1(1) + b.N1(2), b.N2(1) + b.N2(2)], true);
        %
    end
    
    title(cell2mat(Factor_label2(i)));
    
    
    if i ==1
        set(gca,'XTick',[1 2],'XTickLabels',{'W', 'M'});
    elseif i ==2
        set(gca,'XTick',[1 2],'XTickLabels',{'Right', 'Left'});
    elseif i ==3
        set(gca,'XTick',[1:1:4],'XTickLabels',{'Science', 'Social', 'Art', 'Other'});
    elseif i ==4
        set(gca,'XTick', [1:1:7],'XTickLabels',{'0', '<1M', '1M', '>1M', '1W', '>1W', '1'});
    elseif i ==5
        set(gca,'XTick',[1 2 3],'XTickLabels',{'No', 'A little', 'Yes'});
    end
    
    TextLabel =0;
    clear you;
    
    
end
han = axes(fig, 'visible', 'off');
han.YLabel.Visible='on';
ylabel(han,'Proportion');

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig1_Demographic2.jpg')],'-djpeg','-r600');
end

%% For categorical data
% Stat pour données ordinales (un peu, beaucoup... : kruskal_wallis ordinal
% Stat pour données qualitatives :chi square (prop_test)
Init;

Factor_cate = [Data.UsedToEnigmas];
Factor_label2 = {'Used to problem solving'};

figure
decalage = 0;
higherBar = 1;


% Define variables
for i =1:2
    subplot(1,2,i)
    
    Factor = Factor_cate;
    
    if i==1
        you.Wake = Factor(WakeAll); you.N1 = Factor(N1All);
        map = Color2;
        colormap(map);
        [a,b,c] = nhist(you, 'p','Proportion', 'fsize', 18, 'linewidth', 4, 'smooth', 'color','colormap')
        
        groups = [zeros(1,length(b.Wake)-2) ones(1,length(b.N1)-2)];
        data = [b.Wake(1:end-2) b.N1(1:end-2)];
        p(i) = kruskalwallis(data, groups, 'off');
        
    elseif i ==2
        you2.Solvers = Factor(AllSol); you2.NonSolvers = Factor(AllNonSol);
        map = ColorSolver;
        colormap(map);
        [a,b,c] = nhist(you2, 'p','Proportion', 'fsize', 18, 'linewidth', 4, 'smooth', 'color','colormap')
        
        groups = [zeros(1,length(b.Solvers)-2) ones(1,length(b.NonSolvers)-2)];
        data = [b.Solvers(1:end-2) b.NonSolvers(1:end-2)];
        p = kruskalwallis(data, groups, 'off');
        
    end
    ylabel('Proportion');
    
    set(gca,'XTick',[1 2 3],'XTickLabels',{'No', 'A little', 'Yes'});
    
    
end

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Enigma.jpg')],'-djpeg','-r600');
end
