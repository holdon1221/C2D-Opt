function Fig2B_i

clear all

% Input time-series data for follicular phase hormones

% E2 data: time in hours and corresponding hormone levels
E2timeF = [13.5, 15.6, 17.5, 19.4, 21.5, 23.5, 25.4, 27.6, 29.4, 31.5, 33.5,...
                  35.5, 37.5, 39.5, 41.4, 43.4, 45.3, 47.4, 49.3, 51.6, 53.5, 55.3];  
E2dataF = [72.8,76.1,79.4,78.8,74.9,80.4,78.2,73.2,75.3,82.5,82.5,84.1, ...
                  88.5,76.8,81.8,86.7,101,90.5,98.8,92.1,87.1,88.7];
% P4 data
P4timeF = [13.6, 15.6, 17.6, 19.6, 21.3, 23.6, 25.5, 27.4, 29.4, 31.4,...
           33.4, 35.3, 37.5, 39.4, 41.2, 43.2, 45.3, 47.3, 49.4, 51.3, 53.3, 55.4];   
P4dataF = [2.42,2.16,2.17,1.86,2.02,1.79,2.01,2.09,2.46,2.43,2.46,...
           2.35,2.29,2.22,2.18,2.03,2.03,1.97,2.11,2.40,2.53,2.43];

% LH data
LHtimeF = [13.6, 15.7, 17.6, 19.6, 21.6, 23.6, 25.5, 27.6, 29.5,...
           31.5, 33.6, 35.6, 37.5, 39.5, 41.6, 43.6, 45.6, 47.5, 49.6, 51.6, 53.5, 55.5];
LHdataF = [5.10,5.85,5.42,5.64,4.32,4.32,4.10,4.37,4.94,4.55,...
           5.17,5.21,5.74,5.30,5.09,5.13,5.66,5.40,4.96,5.27,4.83,4.70];

% FSH data
FSHtimeF = [13.5, 15.5, 17.4, 19.5, 21.4, 23.5, 25.6, 27.5, 29.5,...
            31.6, 33.6, 35.4, 37.5, 39.5, 41.4, 43.5, 45.4, 47.5, 49.6, 51.5, 53.6, 55.6];   
FSHdataF = [6.29,6.43,6.46,6.80,6.15,6.15,6.01,5.81,6.25,5.84,...
            5.77,5.90,5.97,6.10,6.07,6.27,5.76,5.72,5.42,5.55,5.18,5.31];

% Select time window of interest 
i = 10;
j = 22;

% Compute mean hormone levels over selected window (used for normalization)
meanE2 = mean(E2dataF(10:22));
meanP4 = mean(P4dataF(10:22));
meanLH = mean(LHdataF(10:22));
meanFSH = mean(FSHdataF(10:22));

% Define parameters for cosine fitting of circadian-like hormone patterns ---
    a1 =       0.2739/meanLH;   % LH parameter for amplitude
    c1 =       18.27;           % LH parameter for phase (in hours)
    d1 =       1;               % LH data normalized mean
    a2 =       0.3797/meanFSH;  % FSH parameter for amplitude
    c2 =       16.46;           % FSH parameter for phase (in hours)
    d2 =       1;               % FSH data normalized mean
    a3 =       7.045/meanE2;    % E2 parameter for amplitude
    c3 =      0.5845;           % E2 parameter for phase (in hours)               
    d3 =       1;               % E2 data normalized mean
    a4 =      0.2303/meanP4;    % P4 parameter for amplitude
    c4 =       8.46;            % P4 parameter for phase (in hours)              
    d4 =       1;               % P4 data normalized mean

 
    xdata = [0:0.1:60]';

    % Compute cosine fits for each hormone
    cosineLH = a1*cos((pi/12)*(xdata - c1)) + d1;
    cosineFSH = a2*cos((pi/12)*(xdata - c2)) + d2;
    cosineE2 = a3*cos((pi/12)*(xdata - c3)) + d3;
    cosineP4 = a4*cos((pi/12)*(xdata - c4)) + d4;
    

    % Plot normalized FSH data with fitted cosine curve
    figure(1)
    set(gcf, 'Color', 'white')
    p = plot(xdata, cosineFSH, 'k', 'LineWidth',1);
    hold on
    p = plot(FSHtimeF(i:j), FSHdataF(i:j)./meanFSH, 'ro', 'MarkerSize', 5,'linewidth', 3);
    ylabel('Norm $FSH$', 'Interpreter', 'latex');
    set(gca, 'Fontsize', 20)
    xlim([31, 55])
    xticks([31 55]);
    xtickangle(90);
    xticklabels({'7:00',  '7:00'});
    xlabel('t [hours]', 'Interpreter', 'latex');
    ylim([0.8, 1.2])
    yticks([0.8 1 1.2]);

    % Plot normalized LH data with fitted cosine curve
    figure(2)
    set(gcf, 'Color', 'white')
    p = plot(xdata, cosineLH, 'k', 'LineWidth',1);
    hold on
    p = plot(LHtimeF(i:j), LHdataF(i:j)./meanLH, 'ro', 'MarkerSize', 5,'linewidth', 3);
    ylabel('Norm $LH$', 'Interpreter', 'latex');
    set(gca, 'Fontsize', 20)
    xlim([31, 55])
    xticks([31 55]);
    xtickangle(90);
    xticklabels({'7:00',  '7:00'});
    xlabel('t [hours]', 'Interpreter', 'latex');
    ylim([0.8, 1.2])
    yticks([0.8 1 1.2]);

    % Plot normalized E2 data with fitted cosine curve
    figure(3)
    set(gcf, 'Color', 'white')
    p = plot(xdata, cosineE2, 'k', 'LineWidth',1);
    hold on
    p = plot(E2timeF(i:j), E2dataF(i:j)./meanE2, 'ro', 'MarkerSize', 5,'linewidth', 3);
    ylabel('Norm $E_2$', 'Interpreter', 'latex');
    set(gca, 'Fontsize', 20)
    xlim([31, 55])
    xticks([31 55]);
    xtickangle(90);
    xticklabels({'7:00',  '7:00'});
    xlabel('t [hours]', 'Interpreter', 'latex');
    ylim([0.8, 1.2])
    yticks([0.8 1 1.2]);

    % Plot normalized P4 data with fitted cosine curve
    figure(4)
    set(gcf, 'Color', 'white')
    p = plot(xdata, cosineP4, 'k', 'LineWidth',1);
    hold on
    p = plot(P4timeF(i:j), P4dataF(i:j)./meanP4, 'ro', 'MarkerSize', 5,'linewidth', 3);
    ylabel('Norm $P_4$', 'Interpreter', 'latex');
    set(gca, 'Fontsize', 20)
    xlim([31, 55])
    xticks([31 55]);
    xtickangle(90);
    xticklabels({'7:00',  '7:00'});
    xlabel('t [hours]', 'Interpreter', 'latex');
    ylim([0.8, 1.2])
    yticks([0.8 1 1.2]);

end

