% Plot the trial number of the Insight (for post solvers)


Init;
decalage = 0; higherBar = 0;

figure

x=[0.25+zeros(length(AhaMoment_Wake), 1); 0.75+zeros(length(AhaMoment_N1),1); 1.25+zeros(length(AhaMoment_N2),1)]
y = [AhaMoment_Wake; AhaMoment_N1; AhaMoment_N2];
h = beeswarm(x, y, 'sort_style', 'up', 'overlay_style', 'ci', 'dot_size', 3, 'colormap', Color3N2, 'MarkerFacealpha', 0.5,'MarkerEdgeColor','k');
Labels = 'Trials';
TextLabel = 'Eureka: Trial Number';
p= anova1(padcat(AhaMoment_Wake, AhaMoment_N1, AhaMoment_N2), [],'off');
sigstar([0.25, 1.25], p ,[], higherBar)
% title ('Solvers \it{Post}');

Design;

Labels = Labels3N2;
set(gca,'XTick',[0.25 0.75 1.25],'XTickLabels',Labels,'FontWeight', 'bold');

if SaveFig ==1
    print(gcf,[FigPath,filesep, sprintf('Fig_EurekaMoment.jpg')],'-djpeg','-r600');
end
%%
% Scatter time N1 vs tradeoff before aha
figure
timeN1 = Data.timeN1Edison (Data.InsightPost ==1 & Data.HypnoEdisonN2 ==0);
tradeoff= Data.TradeOff_beforeAha (Data.InsightPost ==1& Data.HypnoEdisonN2 ==0);
AhaMoment = Data.AhaMoment(Data.InsightPost ==1 & Data.HypnoEdisonN2 ==0)
scatter(tradeoff, AhaMoment)


% Fig trade off before Eureka Sol vs NonSol: can we find a trace of an
% eureka before the eureka moment?

Init
decalage = 0; higherBar = 0;

figure

x=[0.25+zeros(length(Data.TradeOff_Start(AllSol)), 1); 0.75+zeros(length(Data.TradeOff_Start(AllNonSol)),1)]
y = [Data.TradeOff_Start(AllSol); Data.TradeOff_Start(AllNonSol)];
h = beeswarm(x, y, 'sort_style', 'up', 'overlay_style', 'ci', 'dot_size', 3, 'colormap', ColorSolver, 'MarkerFacealpha', 0.5,'MarkerEdgeColor','k');

Labels = 'Trial';
TextLabel = 'Trade-off before Eureka';
[h,p]= ttest2(Data.TradeOff_Start(AllSol), Data.TradeOff_Start(AllNonSol));
sigstar([0.25 0.75], p ,[], higherBar)

Design;

Labels = LabelsSol;
set(gca,'XTick',[0.25 0.75],'XTickLabels',Labels,'FontWeight', 'bold');


if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_TradeOff_Before_Eureka.jpg')],'-djpeg','-r600');
end