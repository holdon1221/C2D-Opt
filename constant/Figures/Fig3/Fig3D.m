% ============================================================ 
% Fig3D.m
%
% PURPOSE:
%   Plot total optimized EE and DNG doses (constant administration) across 24 dosing times 
%   (1:00–24:00), based on a 28-day treatment cycle (21 days with dosing,
%   followed by 7 days without treatment).
%   Uses dual y-axes to display EE and DNG doses side-by-side per dosing time.
%   This figure illustrates how optimal total dose varies with intake time.
% ============================================================

function Fig3D

% Define the 24-hour dosing times (from 1:00 AM to 12:00 AM)
days = 1:24;  

% Total EE dose per dosing time
% Values derived from optimized simulations; converted from ng to µg
data1 = [630000			% Total EE dose (optimized) for 1:00 dosing time	
630000				% Total EE dose (optimized) for 2:00 dosing time
630000				% Total EE dose (optimized) for 3:00 dosing time
630000				% Total EE dose (optimized) for 4:00 dosing time
630000				% Total EE dose (optimized) for 5:00 dosing time
630000				% Total EE dose (optimized) for 6:00 dosing time
604800				% Total EE dose (optimized) for 7:00 dosing time
604800				% Total EE dose (optimized) for 8:00 dosing time
579600				% Total EE dose (optimized) for 9:00 dosing time
604800				% Total EE dose (optimized) for 10:00 dosing time
592200				% Total EE dose (optimized) for 11:00 dosing time
592200				% Total EE dose (optimized) for 12:00 dosing time
617400				% Total EE dose (optimized) for 13:00 dosing time
592200				% Total EE dose (optimized) for 14:00 dosing time
630000				% Total EE dose (optimized) for 15:00 dosing time
617400				% Total EE dose (optimized) for 16:00 dosing time
630000				% Total EE dose (optimized) for 17:00 dosing time
630000				% Total EE dose (optimized) for 18:00 dosing time
617400				% Total EE dose (optimized) for 19:00 dosing time
630000				% Total EE dose (optimized) for 20:00 dosing time
630000				% Total EE dose (optimized) for 21:00 dosing time
630000				% Total EE dose (optimized) for 22:00 dosing time
630000				% Total EE dose (optimized) for 23:00 dosing time
630000]/1000;			% Total EE dose (optimized) for 24:00 dosing time

% Total DNG dose per dosing time
% Values derived from optimized simulations; converted from ng to µg
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


% Start figure and use dual y-axes to compare EE and DNGfigure;

% Plot EE doses on left y-axis (blue bars)
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

% Plot DNG doses on right y-axis (magenta bars)
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

