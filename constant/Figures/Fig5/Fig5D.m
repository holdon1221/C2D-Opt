% ============================================================
% Fig5D.m
%
% PURPOSE:
%   This script compares total EE and DNG doses (in µg) for optimized constant vs optimized
%   nonconstant administration across 24 dosing times (1:00–24:00), based on a 
%   28-day treatment cycle (21 days with dosing, followed by 7 days without treatment).
%   Constant administration assumes the same daily dose each day, but the dose level
%   itself is optimized. Nonconstant administration uses optimized
%   day-specific doses and possibly fewer intake days.
%   Dual y-axes show EE and DNG dose per dosing time, with overlaid bars allowing
%   direct comparison of the two optimized strategies.
% ============================================================

function Fig5D

% Define 24 dosing times (1:00 - 24:00)
days = 1:24;  

% Total EE dose (optimized) for nonconstant administration based on a 28-day treatment cycle, 
% across dosing times 1:00 to 24:00.
data1 = [259.23, 		% Total EE dose (optimized) for 1:00 dosing time
	259.23, 		% Total EE dose (optimized) for 2:00 dosing time
	209.58, 		% Total EE dose (optimized) for 3:00 dosing time
	209.58, 		% Total EE dose (optimized) for 4:00 dosing time
	209.58, 		% Total EE dose (optimized) for 5:00 dosing time
	209.58, 		% Total EE dose (optimized) for 6:00 dosing time
	209.58, 		% Total EE dose (optimized) for 7:00 dosing time
	209.58, 		% Total EE dose (optimized) for 8:00 dosing time
	193.78, 		% Total EE dose (optimized) for 9:00 dosing time
	193.78, 		% Total EE dose (optimized) for 10:00 dosing time
    193.78, 		% Total EE dose (optimized) for 11:00 dosing time
	193.78, 		% Total EE dose (optimized) for 12:00 dosing time
	193.78, 		% Total EE dose (optimized) for 13:00 dosing time
	210.83, 		% Total EE dose (optimized) for 14:00 dosing time
	210.83, 		% Total EE dose (optimized) for 15:00 dosing time
	210.83, 		% Total EE dose (optimized) for 16:00 dosing time
	214.47, 		% Total EE dose (optimized) for 17:00 dosing time
	229.22, 		% Total EE dose (optimized) for 18:00 dosing time
	229.22, 		% Total EE dose (optimized) for 19:00 dosing time
	240.09, 		% Total EE dose (optimized) for 20:00 dosing time
    240.09, 		% Total EE dose (optimized) for 21:00 dosing time
	259.88, 		% Total EE dose (optimized) for 22:00 dosing time
	259.88, 		% Total EE dose (optimized) for 23:00 dosing time
	259.88];		% Total EE dose (optimized) for 24:00 dosing time

% Total DNG dose (optimized) for nonconstant administration based on a 28-day treatment cycle,
% across dosing times 1:00 to 24:00.
data2 = zeros(1,24);

% Total EE dose (optimized) for constant administration based on a 28-day treatment cycle, 
% across dosing times 1:00 to 24:00.
data3 = [630000, 630000, 630000, 630000, 630000, 630000, 604800, 604800, 579600, 604800, ...
         592200, 592200, 617400, 592200, 630000, 617400, 630000, 630000, 617400, 630000, ...
         630000, 630000, 630000, 630000]/1000;

% Total DNG dose (optimized) for constant administration based on a 28-day treatment cycle, 
% across dosing times 1:00 to 24:00.
data4 = [15120000, 13440000, 11760000, 10080000, 9240000, 8400000, 8400000, 7560000, 7560000, ...
         6720000, 6720000, 6720000, 6720000, 7560000, 7560000, 8400000, 9240000, 10080000, ...
         11760000, 12600000, 13440000, 14280000, 13440000, 12600000]/1000;

% Number of intake days under nonconstant administration, across dosing times 1:00 to 24:00.
data5 = [10, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, ...
         11, 10, 10, 10];

% ====== Plot doses under constant and nonconstant strategies ======
figure(1);
hold on; 

% Plot nonconstant EE (solid blue bars)
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

% Plot nonconstant DNG (solid magenta bars)
yyaxis right;
bar2 = bar(days + 0.2, data2, 0.4, 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');  
ylim([0 16000]);     % 4500
yticks([0 16000]);
yticklabels([0 16000]);
ylabel('DNG [$\mu$g]', 'Interpreter', 'latex');
ax = gca;  
ax.YColor = [1, 0, 1];  

% Plot constant EE (blue transparent overlay)
yyaxis left;
bar3 = bar(days - 0.2, data3, 0.4, 'FaceColor', [0, 0, 1], 'EdgeColor', 'black', 'FaceAlpha', 0.3, 'LineWidth', 0.5);

% Plot constant DNG (magenta transparent overlay)
yyaxis right;
bar4 = bar(days + 0.2, data4, 0.4, 'FaceColor', [1, 0, 1], 'EdgeColor', 'black', 'FaceAlpha', 0.2, 'LineWidth', 0.5);

xlim([0 25]);
xticks([1 6 12 18 24]);
xticklabels([1 6 12 18 24]);
xtickangle(90)
xlabel('Dosing time [hours]', 'Interpreter', 'latex');

% ====== Plot number of intake days (nonconstant administration only) ======
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
