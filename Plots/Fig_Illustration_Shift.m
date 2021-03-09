Init;
decalage = 1; higherBar = 0;

fig = figure;
subplot(1,2,1)

load('C:\Users\Célia\Desktop\WagnerEdison project\NRT_07_07_2020_Scripts\Results\10OL_27-Jul-2020\NRT_10OL_post.mat');
Post_T= struct2table(Post);Post_T.Time = smoothdata(Post_T.Time);

bar(Post_T.Time(1:270), 'FaceColor', verteau);
title('Solver', 'color', verteau);
 TextLabel = 'Solving Time (secs)';
Labels = 'Trials'; 
ahamoment = findchangepts(Post_T.Time)
xline(ahamoment, '--', 'color','k', 'linewidth',3)

% xline(271, '--', 'color', 'k', 'linewidth', 3)
Design;
xlabel('Trial Number', 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color', Basiccolor);
%%

%Solver

subplot(1,2,2)


load('C:\Users\Célia\Desktop\WagnerEdison project\NRT_07_07_2020_Scripts\Results\01HC_20-Jul-2020\NRT_01HC_post.mat');
Post_T= struct2table(Post); Post_T.Time = smoothdata(Post_T.Time);

bar(Post_T.Time(1:270), 'FaceColor', red);
title('Non Solver', 'color', red)
Labels = 'Trials'; TextLabel = 'Solving Time (secs)';
ahamoment = findchangepts(Post_T.Time)
xline(ahamoment, '--', 'color','k', 'linewidth',3)
% xline(271, '--', 'color', 'k', 'linewidth', 3)
Design;
xlabel('Trial Number', 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color', Basiccolor);

% Give common xlabel, ylabel and title to your figure
% han=axes(fig,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Solving Time (secs)');
% xlabel(han,'Trial Number');

% h = legend('Insight', 'Implicite'); 

% legend boxoff

Design;

if SaveFig ==1
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Illustration_shift.jpg')],'-djpeg','-r600');
end
