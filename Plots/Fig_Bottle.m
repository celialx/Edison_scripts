% Pense a rajouter les sujets bug en Post mais qui ont fait la partie bouteille : 13AB & 66LY

% Fig: Nb total repet bottle
% Init;
% 
% figure
% decalage = 0;
% c=cbrewer('div', 'Spectral', 3);
% bar(1, Nb_Ball_Single, 'FaceColor', c(1,:), 'LineWidth', 2, 'EdgeColor', Basiccolor); hold on; bar(2, Nb_Ball, 'FaceColor', c(2,:), 'LineWidth', 2, 'EdgeColor', Basiccolor);
% 
% str1 = sprintf('N = %d',Nb_Ball_Single); Text1 =text(1, max(Nb_Ball_Single+3), str1, 'color', Basiccolor, 'FontWeight', 'bold');
% str2 = sprintf('N = %d',Nb_Ball); Text2 =text(2, max(Nb_Ball+3), str2, 'color', 0.1*c(2,:) ,'FontWeight', 'bold');
% 
% TextLabel = 'Number of bottle dropped';
% Labels = {'Once','Several'};
% 
% set(Text1, 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color',Basiccolor, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
% set(Text2, 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color',Basiccolor, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
% 
% Design;
% 
% if SaveFig ==1
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
%     print(gcf,[FigPath,filesep, sprintf('Fig_General_Bottle.jpg')],'-djpeg','-r600');
% end

% Fig : A: timeN1 before drop percentage, B: timeN1before drop, C: meta
% sleep


Init;
decalage = 1; higherBar = 1;
c=cbrewer('div', 'Spectral', 3);

subplot(2,2,1) % 1 prevent from going into N2
% drop_beforeN2 = (length(stadeBall_All)-sum(stadeBall_All==2))/length(stadeBall_All)*100
% bar(1, drop_beforeN2,'FaceColor', [0 0.75 0.75], 'LineWidth', 2, 'EdgeColor', Basiccolor);
% str1 = sprintf('%d %',round(drop_beforeN2)); Text1 =text(1, drop_beforeN2+5, str1, 'color', Basiccolor, 'FontWeight', 'bold');
% hold on
% TextLabel = 'Bottle drop (%)';
% Labels = {'Before N2'};
% Design; 

subplot(2,2,1) % % fell N1 vs wake

bar(1,nnz(timeN1bfBall)/length(timeN1bfBall)*100,'FaceColor', N1color, 'LineWidth', 2, 'EdgeColor', N1edge)
str1 = sprintf('%d %',nnz(timeN1bfBall)/length(timeN1bfBall)*100); Text1 =text(1, max(nnz(timeN1bfBall)/length(timeN1bfBall)*100+2), str1, 'color', Basiccolor, 'FontWeight', 'bold');
hold on

bar(2,(length(timeN1bfBall)-nnz(timeN1bfBall))/length(timeN1bfBall)*100,'FaceColor', wakecolor, 'LineWidth', 2, 'EdgeColor', wakeedge)
str2 = sprintf('%d %',(length(timeN1bfBall)-nnz(timeN1bfBall))/length(timeN1bfBall)*100); Text2 =text(2, (length(timeN1bfBall)-nnz(timeN1bfBall))/length(timeN1bfBall)*100+2, str2, 'color', Basiccolor, 'FontWeight', 'bold');
hold on

set(Text1, 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color',Edge2(2,:), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
set(Text2, 'FontSize', 20, 'FontWeight','bold','FontName','Helvetica','Color',Edge2(1,:), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');


TextLabel = 'Bottle drop after N1 (%)';
Labels = {'Yes'; 'No'};

[h,p, chi] = prop_test([nnz(timeN1bfBall), length(timeN1bfBall)-nnz(timeN1bfBall)], [length(timeN1bfBall),length(timeN1bfBall)], true);
higherBar = 1; % Pour passer au dessus des surtitres
sigstar([1,2],p, [], higherBar);
Design;

subplot(2,2,2)
gyt_PiratePlot(1,timeN1bfBall/60,0.5,2,'y',N1color,'y');
hold on
Labels = {'Micro sleep'};

TextLabel = 'Time of N1 before drop (min)';
set(gca, 'xlim', [0.3 1.5]);
Design;

% subplot(2,3,4)
% 
% load([path_local 'Tf_around_drops.mat'])
% freqs=TFRhann.freq;
% times=TFRhann.time;
% temp_toplot=squeeze(mean(all_TF_drops(:,2,:,:),1));
% h=simpleTFplot(temp_toplot,freqs,times,0,0);
% caxis([-4 3])
% colorbar;
% 
% xlabel('Time from Drop (s)')
% ylabel('Frequency (Hz)')
% title(sprintf('N=%g drops',size(all_TF_drops,1)))
% format_fig;
% clear xlim

subplot(2,2,3)
save_path='C:\Users\Célia\Desktop\WagnerEdison project\Analyses_Results\LocalSleep\';
load([save_path 'Tf_around_drops.mat'])
freqs=TFRhann.freq;
times=TFRhann.time;
FOI=[1 4];
COI=2;
 temp_toplot=squeeze(mean(all_TF_drops(:,COI,freqs>=FOI(1) & freqs<=FOI(2),:),3));
simpleTplot(times,temp_toplot,0,'k',0,'-',0.5,1,0,1,1);
clear xlim;
xlim([-60 30])
format_fig;
xlabel('Time from Drop (s)')
ylabel('Power (dB)')
title(sprintf('N=%g drops',size(all_TF_drops,1)))

subplot(2,2,4)

Labels = {'Yes'; 'No'};
TextLabel = 'Hypnagogia (%)';

h = bar(HypnaBottle/Nb_Bottle*100, 'FaceColor', [0.25 0.25 0.25])
hold on
bar(2, HypnaGeneral/Nb_NoBottle*100, 'FaceColor', [0.5 0.5 0.5])
hold on
Design;

[h,p,chi] = prop_test([HypnaBottle HypnaGeneral], [Nb_Bottle Nb_NoBottle], 'true');
higherBar = 0;
sigstar([1,2],p, [], higherBar);
xlabel('Bottle dropped');


% 
% subplot(1,3,3)
% MetaSleep = nnz(Data.MetaSleep(Data.Fall>0 & Data.HypnoEdisonN1 ==0));
% Fall_Nb = nnz(Data.Fall>0 & Data.HypnoEdisonN1 ==0);
% MetaSleep_Percent = (MetaSleep/Fall_Nb)*100;
% 
% bar(1,MetaSleep_Percent, 'FaceColor',Color2(1,:), 'EdgeColor', Edge2(1,:), 'LineWidth', 2);
% 
% TextLabel = 'Meta-sleep at bottle drop (%)';
% Labels = {'Wake'};
% yticks([0:20:100]);

Design;

if SaveFig ==1
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.70 1]);
    print(gcf,[FigPath,filesep, sprintf('Fig_Bottle.jpg')],'-djpeg','-r600');
end