%% PVT Results
 Init;

% Pense à enlever le premier PVT pour le mean
close all

% Fig PVT
figure
decalage = 1;
higherBar = 0;

Factor = {Wake N1 N2};
Var = {Data.meanPVT_Pre, Data.meanPVT_Post}


for i = 1:6
    
    if i ==1
        subplot(1,2,1)
        title('\it{Pre}')
        bplot(Var{i}(Factor{i}),1, 'color',Color3N2(1,:));
        
    elseif i ==2
        bplot(Var{i-1}(Factor{i}),2, 'color',Color3N2(2,:));
    elseif i ==3
        bplot(Var{i-2}(Factor{i}),3, 'color',Color3N2(3,:));
        
    elseif i ==4
        subplot(1,2,2);
        title('\it{Post}');
        bplot(Var{i-2}(Factor{i-3}),1, 'color',Color3N2(1,:));
        xlim_curr = get(gca,'xlim');
        
    elseif i ==5
        bplot(Var{i-3}(Factor{i-3}),2, 'color',Color3N2(2,:));
    elseif i ==6
        bplot(Var{i-4}(Factor{i-3}),3, 'color',Color3N2(3,:));
        
    end
    
    
    TextLabel = 'Reaction Time (ms)'
    Labels = Labels3N2; Design;
    hold on;
    if i ==3
        x = padcat(Data.meanPVT_Pre(Wake),Data.meanPVT_Pre(N1),Data.meanPVT_Pre(N2));
        p = anova1(x, [], 'off');
        sigstar([1 3],p, [], higherBar);
        
    elseif i == 6
        x = padcat(Data.meanPVT_Post(Wake),Data.meanPVT_Post(N1),Data.meanPVT_Post(N2));
        [p, ANOVATAB, stats] = anova1(x, [], 'off');
        a=multcompare(stats, 'display','off')
        
        sigstar([1 2],a(1,6), [], higherBar);
        sigstar([1 3],a(2,6), [], higherBar);
        sigstar([2 3],a(3,6), [], higherBar);
        
    end
    
end


if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_PVT.jpg')],'-djpeg','-r600');
end