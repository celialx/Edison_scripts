function h = do_linepersuj(data) 
%%% create a single subjects line plot graphic per condition
%%% laure bottemanne
%%% last update 18dec2016

sujNb = size(data.data,1);
x = 1:size(data.data,2);
abs = [x(1)-0.5 x(end)+0.5];

h = figure('Color',[1 1 1]);

axes1 = axes('Parent',h);

% plot1 = 'plot1 = plot(' ; %'plot1 = plot(';
for s = 1:sujNb
    if s == sujNb
        plot1 = [plot1 ' data.x, data.data(' num2str(s) ',:) '];
    else % s < sujNb
        plot1 = [plot1 ' data.x, data.data(' num2str(s) ',:), '];
    end; % if s == sujNb
end; % for s
% plot1 = [plot1 ');'];
% eval(plot_str)

%%% figure display
set(plot1, 'MarkerSize', 10, 'LineWidth',2, 'Marker', 'x', 'LineStyle','--'); 
set(axes1, 'xlim', abs);
set(axes1,...
    'box', 'off', ...
    'YGrid','off',...
    'XGrid', 'off', ...
    'XTickLabel',{char(data.conditions)},...
    'XTick',x,... 
    'FontSize', 20,...
    'FontName','Gill Sans MT');
xlabel(data.xfactor);
ylabel(data.var);