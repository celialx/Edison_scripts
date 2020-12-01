% Training and pre perf per Group (Wake vs N1)

Init;

decalage =1;
higherBar=0;

for f =1:2
    figure
    if f ==1
        sgtitle('Training', 'FontWeight', 'bold', 'FontSize', 20, 'FontName', 'Helvetica')
        Var = {Data.Correct_Practice(Wake), Data.Correct_Practice(N1), Data.Correct_Practice(N2), Data.RT_Practice(Wake), Data.RT_Practice(N1), Data.RT_Practice(N2)};
    elseif f ==2
        sgtitle ('\it{Pre}', 'FontWeight', 'bold', 'FontSize', 20, 'FontName', 'Helvetica')
        Var = {Data.Correct_Pre(Wake), Data.Correct_Pre(N1), Data.Correct_Pre(N2), Data.RT_Pre(Wake), Data.RT_Pre(N1), Data.RT_Pre(N2)};
    end
    
    
    for i = 1:length(Var)
        Labels = Labels3N2;
        
        if i ==1 || i == 2 || i == 3
            subplot(1,2,1)
            bplot(Var{i},i,'color', Color3N2(i,:), 'LineWidth', 2.5);
            
            TextLabel = 'Accuracy (%)'; Design;
            if i ==3
                [p, stats] = anova1(padcat(Var{1}, Var{2}, Var{3}), [], 'off');
                sigstar([1,3], p, [], higherBar);
            end
        elseif i == 4 || i == 5 || i ==6
            subplot(1,2,2)
            bplot(Var{i},i-2,'color', Color3N2(i-3,:), 'LineWidth', 2.5);
            
            TextLabel = 'Solving Time (secs)';Design;
            
            if i ==6
                [p, stats] = anova1(padcat(Var{4}, Var{5}, Var{6}), [], 'off');
                sigstar([2,4], p, [], higherBar);
            end
        end
        
        hold on;
        Design;
        
    end
    
    
end

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_InitialPerf.jpg')],'-djpeg','-r600');
end
