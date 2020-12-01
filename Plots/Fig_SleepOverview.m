%-----------------------------------------------------------------
% Figure - Overview sleep statistics
%------------------------------------------------------------------

close all;
Init; % Contain all the variables of interest

figure;
decalage = 1; % for subplot, shift y label to the left (away from the axis)

% Number of subjects in each vigilance state
subplot(1,2,1);
Var = {Nb_WakeAll, Nb_N1All,Nb_N2All};

for i =1:3
    
    bar(i, Var{i}, 'FaceColor',Color3N2(i,:), 'EdgeColor', Edge3N2(i,:), 'LineWidth', 2);
    text(i, max(Var{i})+1.5, num2str(Var{i}),'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color',Edge3N2(i,:))
    hold on
end

TextLabel = 'Number of Subjects';
Labels = Labels3N2; Design;
text(3.25,ylim(2),['N = ' num2str(Nb_Sub)],'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color',[0.2 0.2 0.2] )
Design;


% Time spent in each stage

subplot(1,2,2)

timeN1_AASM =(Data.NbEpochsN1Edison(N1All)*30)/60; % Conventional scoring
timeN1_Micro_N1=(Data.timeN1Edison(N1All))/60; % BERN scoring
timeN1_Micro_Wake = (Data.timeN1Edison(WakeAll))/60;

gyt_PiratePlot(1,timeN1_AASM,0.3,2,'y',N1color,'y') % violin plot
hold on
gyt_PiratePlot(2,timeN1_Micro_N1,0.3,2,'y',N2color,'y') % violin plot
hold on
gyt_PiratePlot(3,timeN1_Micro_Wake,0.3,2,'y',wakecolor,'y') % violin plot

Labels = {'AASM', 'BERN', 'BERN\_Wake'} ;
TextLabel='Time spent in N1 (min)';
set(gca, 'xlim', [0.5 3.5]);
Design;
box off

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig1_SleepOverview.jpg')],'-djpeg','-r600');
end
