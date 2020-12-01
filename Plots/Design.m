
% General_Design
Basiccolor =  [0.2 0.2 0.2];

if strcmp(Labels, 'Trials')  
else
set(gca,'XTick',[1:length(Labels)],'XTickLabels',Labels,'FontWeight', 'bold');
end

if TextLabel ==0
    ylab = [];
else
ylab = ylabel(TextLabel, 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color', Basiccolor);
end

if exist('ylab') & decalage ==1
set(ylab, 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
elseif decalage ==0
set(ylab, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);

end
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');

% 
% set(gca,'color','none') %transparent background

set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 3, ...
'DefaultAxesXColor', Basiccolor, ...
'DefaultAxesYColor', Basiccolor, ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 20, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 20, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

box off;

% % Ouvre la figure en plein écran
% fig=gcf;
% fig.Units='normalized';
% fig.OuterPosition=[0 0 1 1];

