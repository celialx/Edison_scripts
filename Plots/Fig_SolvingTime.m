% Fig on Solving Time, accuracy, trade off before and after Eureka

Init;

% Fig1 Solving Time before and after Eureka for sol vs non sol
% Fig2: Speed-Accuracy Tradeoff for sol vs non sol
% Figure 3: Tradeoff start of post session vs tradeoff 20 trials before Eureka

decalage =1; higherBar = 0;
cte = 0;

Var1 = {Data.RT_beforeAha(AllSol), Data.RT_fromAha(AllSol), Data.RT_beforeAha(AllNonSol), Data.RT_fromAha(AllNonSol)}; % Figure 1: Solving Time
Var2 = {Data.TradeOff_beforeAha(AllSol), Data.TradeOff_fromAha(AllSol), Data.TradeOff_beforeAha(AllNonSol), Data.TradeOff_fromAha(AllNonSol)}; % Figure 2: TradeOff
% Var3 = {Sol_Start, Sol_beforeAha, Non_Start, NonSol_beforeAha}; % Figure 3: Tradeoff start of post session vs tradeoff 20 trials before Eureka

for fig =1:2
    
    figure
    if fig ==1;
        Var = Var1;
    elseif fig ==2
        cte = 0;
        Var = Var2;
    elseif fig ==3
        cte = 0;
        Var = Var3;
    end
    
    for i =1:2
        
        if i ==2
            cte = cte+1;
        end
        
        subplot(1,2,i);
        
        N = length(Var{i+cte}); X = linspace(0,pi*3,1000); C = linspecer(N);
        
        for ii=1:N
            plot([Var{i+cte}(ii), Var{i+cte+1}(ii)]','LineWidth', 2, 'color', C(ii,:), 'linestyle', '--')
            hold on;
        end
        
        if fig ==1
            TextLabel = 'Solving Time (secs)';
            
        elseif fig ==2
            TextLabel = 'Speed-Accuracy Trade-off';
        end
        
        if fig ==3
            Labels = {'Start', 'Before_Aha'};
            
        else
            Labels = {'Before', 'After'};
        end
        Design;
        
        
        [h,p] = ttest(Var{i+cte}, Var{i+cte+1});
        sigstar([1:2], p, [], higherBar);
        
        set(gca, 'xlim', [0.5 2.5]);
        if i ==1
            title('Solvers', 'color', verteau);
        elseif i ==2
            title('Non Solvers', 'color', red);
        end
    end
    
    if fig ==1  && SaveFig ==1
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
        print(gcf,[FigPath,filesep, sprintf('Fig_SolvingTime.jpg')],'-djpeg','-r600');
    elseif fig ==2 && SaveFig ==1
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
        print(gcf,[FigPath,filesep, sprintf('Fig_TradeOff.jpg')],'-djpeg','-r600');
    end
end


%% Percentage correct before and after Eureka
Init;
figure;

decalage = 1; higherBar = 0;
Var = {Data.Correct_beforeAha(AllSol), Data.Correct_fromAha(AllSol), Data.Correct_beforeAha(AllNonSol),Data.Correct_fromAha(AllNonSol)};
Color = [verteau2; verteauedge2; red; rededge];
cte = 0;

for i =1:2
    
    subplot(1,2,i);
    
    if i ==2
        cte = cte+1;
    end
    
    gyt_PiratePlot(1,Var{i+cte},0.35,0.75,'y',Color(i+cte,:),'n')
    hold on
    gyt_PiratePlot(2,Var{i+cte+1},0.35,0.75,'y',Color(i+cte+1,:),'n')
    
    TextLabel = 'Accuracy (%)';
    Labels = {'Before', 'After'};
    
    if i ==1
        title('Solvers', 'color', verteau);
    elseif i ==2
        title('Non Solvers', 'color', red);
    end
    
    
    [h,p] = ttest(Var{i+cte}, Var{i+cte+1});
    sigstar([1:2], p, [], higherBar);
    
    set(gca, 'xlim', [0.5 2.5]);
    set(gca, 'ylim', [40 110]);
    yticks([40:10:100]);
    Design;
    
end

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Accuracy.jpg')],'-djpeg','-r600');
end

%% Figure controle implicite for solvers with no clear drop in solving time

Init
figure
higherBar =1; decalage = 0;
ID = {'20JJ'; '83MB'; '84RD'; '90ED'};

current =0;
for i = 1:length(Data.Code)
    SOI = nnz(strcmp(Data.Code(i), ID));
    
    if SOI>0
        current =current +1;
        Impli_Correct_Rule(current) = Data.Correct_Implicit_Rule(i);
        Impli_Correct_Break(current) = Data.Correct_Implicit_BreakRule(i);
    else
        continue
    end
end

N = length(Impli_Correct_Rule); X = linspace(0,pi*3,1000); C = linspecer(N);

hold off;
var = [Impli_Correct_Rule, Impli_Correct_Break]'
for ii=1:N
    plot([Impli_Correct_Rule(ii), Impli_Correct_Break(ii)]','LineWidth', 2, 'color', C(ii,:), 'linestyle', '--')
    hold on;
end
TextLabel = 'Correct Implicit';
Labels = {'Rule', 'No Rule'};
set(gca, 'xlim', [0.5 2.5]);

Design;

diff = minus(Impli_Correct_Rule, Impli_Correct_Break)
Meandiff = mean(diff);
