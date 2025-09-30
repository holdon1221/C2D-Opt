function Fig3D

days = 1:24;  

data1 = [630000
630000
630000
630000
630000
630000
604800
604800
579600
604800
592200
592200
617400
592200
630000
617400
630000
630000
617400
630000
630000
630000
630000
630000]/1000;

data2 = [15120000
13440000
11760000
10080000
9240000
8400000
8400000
7560000
7560000
6720000
6720000
6720000
6720000
7560000
7560000
8400000
9240000
10080000
11760000
12600000
13440000
14280000
13440000
12600000]/1000;


figure;
yyaxis left;  
bar1 = bar(days - 0.2, data1, 0.4, 'FaceColor', [0, 0, 1]);  
set(gca,'FontSize',18)
ylim([0 700]);
yticks([0 700]);
yticklabels([0 700]);
ylabel('Total EE dose [$\mu$g]', 'Interpreter', 'latex');
ax = gca;  
ax.YColor = [0, 0, 1];  
hold on;

yyaxis right;  
bar2 = bar(days + 0.2, data2, 0.4, 'FaceColor', [1, 0, 1]);  
set(gca,'FontSize',18)
ylim([0 16000]);
yticks([0 16000]);
yticklabels([0 16000]);
ylabel('Total DNG dose [$\mu$g]','Interpreter','latex')
ax.YColor = [1, 0, 1];  

xlim([0 25]);
xticks([1 6 12 18 24]);
xticklabels([1 6 12 18 24]);
xtickangle(90)
xlabel('Dosing time [hours]','Interpreter','latex')

