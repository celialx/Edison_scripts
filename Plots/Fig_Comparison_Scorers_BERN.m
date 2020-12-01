Init;

decalage =1;
higherBar = 0;

load('C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Two_scorers\Edison_Tables_CL\Edison_EEG_CL.mat');
load('C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\Two_scorers\Edison_Tables_DO\Edison_EEG_DO.mat');

Data_CL = Edison_EEG_CL; Data_DO = Edison_EEG_DO;
timeN1_Micro_CL=Data_CL.timeN1Edison;

timeN1_Micro_DO=Data_DO.timeN1Edison;

subplot(1,2,1)
N = length(timeN1_Micro_CL);
X = linspace(0,pi*3,1000);
C = linspecer(N);

hold off;
var = [timeN1_Micro_CL, timeN1_Micro_DO]';

for ii=1:N
    plot([timeN1_Micro_CL(ii), timeN1_Micro_DO(ii)]','LineWidth', 2, 'color', C(ii,:), 'linestyle', '--');
    hold on;
end
TextLabel = 'BERN: Time of N1 (secs)';
Labels = {'CL', 'DO'};

[h,p] = ttest2(timeN1_Micro_CL, timeN1_Micro_DO);
sigstar([1:2], p, [], higherBar);

set(gca, 'xlim', [0.5 2.5]);
Design;


subplot(1,2,2)

gyt_PiratePlot(1,timeN1_Micro_CL,0.3,0.5,'y',N1color,'y');
hold on;
gyt_PiratePlot(2,timeN1_Micro_DO,0.3,0.5,'y',N1color,'y');

Labels = {'CL', 'DO'} ;

TextLabel='BERN: Time of N1 (secs)';
set(gca, 'xlim', [0.5 3.5]);
Design;
box off;

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('BERN_inter_scorers.jpg')],'-djpeg','-r600');
end