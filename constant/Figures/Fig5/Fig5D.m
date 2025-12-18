% This script generates two figures for Fig 5D:
% - Figure 1: Stacked bar plot comparing EE and DNG total doses per time-of-day (1:00-24:00) for constant and nonconstant drug administration.
% - Figure 2: Total number of dosing days for nonconstant administration.

function Fig5D

days = 1:24;  

% Optimal total EE dose for anovulation for each dosing time (1:00–24:00) under nonconstant administration.
data1 = [259.23, 259.23, 209.58, 209.58, 209.58, 209.58, 209.58, 209.58, 193.78, 193.78, ...
         193.78, 193.78, 193.78, 210.83, 210.83, 210.83, 214.47, 229.22, 229.22, 240.09, ...
         240.09, 259.88, 259.88, 259.88];
% Optimal total DNG dose for anovulation for each dosing time (1:00–24:00) under nonconstant administration.
data2 = zeros(1,24);

% Optimal total EE dose for anovulation for each dosing time (1:00–24:00) under constant administration.
data3 = [630000, 630000, 630000, 630000, 630000, 630000, 604800, 604800, 579600, 604800, ...
         592200, 592200, 617400, 592200, 630000, 617400, 630000, 630000, 617400, 630000, ...
         630000, 630000, 630000, 630000]/1000;
% Optimal total DNG dose for anovulation for each dosing time (1:00–24:00) under constant administration.
data4 = [15120000, 13440000, 11760000, 10080000, 9240000, 8400000, 8400000, 7560000, 7560000, ...
         6720000, 6720000, 6720000, 6720000, 7560000, 7560000, 8400000, 9240000, 10080000, ...
         11760000, 12600000, 13440000, 14280000, 13440000, 12600000]/1000;
% Total dosing days for each dosing time (1:00–24:00) under nonconstant administration.
data5 = [10, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, ...
         11, 10, 10, 10];

figure(1);
hold on; 
yyaxis left;
bar1 = bar(days - 0.2, data1, 0.4, 'FaceColor', [0, 0, 1], 'EdgeColor', 'none');  
ylabel('$EE$ [$\mu$g]', 'Interpreter', 'latex');
set(gca,'FontSize',18)
ylim([0 700]);
yticks([0 250 700]);
yticklabels([0 250 700]);
ylabel('EE [$\mu$g]', 'Interpreter', 'latex');
ax = gca;  
ax.YColor = [0, 0, 1];  

yyaxis right;
bar2 = bar(days + 0.2, data2, 0.4, 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');  
ylim([0 16000]);     % 4500
yticks([0 16000]);
yticklabels([0 16000]);
ylabel('DNG [$\mu$g]', 'Interpreter', 'latex');
ax = gca;  
ax.YColor = [1, 0, 1];  

yyaxis left;
bar3 = bar(days - 0.2, data3, 0.4, 'FaceColor', [0, 0, 1], 'EdgeColor', 'black', 'FaceAlpha', 0.3, 'LineWidth', 0.5);

yyaxis right;
bar4 = bar(days + 0.2, data4, 0.4, 'FaceColor', [1, 0, 1], 'EdgeColor', 'black', 'FaceAlpha', 0.2, 'LineWidth', 0.5);

xlim([0 25]);
xticks([1 6 12 18 24]);
xticklabels([1 6 12 18 24]);
xtickangle(90)
xlabel('Dosing time [hours]', 'Interpreter', 'latex');

figure(2)
bar5 = bar(days, data5, 0.5, 'FaceColor', [0.8, 0.8, 0.8]); 
set(gca,'FontSize',18)
ylim([5 10]);
yticks([5 6 7 8 9 10]);
yticklabels([5 6 7 8 9 10]);
ylabel('Total dosing days', 'Interpreter', 'latex');
xlim([0 25]);
xticks([1 6 12 18 24]);
xticklabels([1 6 12 18 24]);
xtickangle(90)
xlabel('Dosing time [hours]', 'Interpreter', 'latex');
