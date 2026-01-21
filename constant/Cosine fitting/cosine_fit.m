% =========================================================================
% cosine_fit.m
%
% PURPOSE (paper correspondence):
%   Fit a simple 24-hour cosine function to *normalized* circadian hormone data,
%   following the workflow described in the paper (Fig. 2B and Methods).
%
%   In the paper, the circadian component is obtained by:
%     1) extracting circadian hormone measurements,
%     2) normalizing by the mean so the average equals 1,
%     3) fitting a cosine curve of the form  a*cos(2*pi*(t - b)) + 1,
%        where:
%          - a = amplitude,
%          - b = acrophase (timing of peak),
%          - t = time (days).
%
%   This script implements the same idea, but uses time in *hours* (x) rather
%   than days (t):
%     cos((pi/12)*(x - c1))  ==  cos(2*pi*(x/24 - c1/24))
%   so c1 is an acrophase expressed in hours (0–24).
%
% OUTPUT:
%   MATLAB 'fit' object f1 with parameters:
%     a1 : fitted amplitude (unitless, because y is normalized)
%     c1 : fitted acrophase (hours)
%     d1 : fitted offset (should be near 1 after normalization)
% =========================================================================
% This script performs a cosinor fit to circadian hormone data (E2, P4, LH, FSH).
% It normalizes the data and fits a cosine curve to capture the rhythmic pattern.

% -------------------------------------------------------------------------
% Main entry point
% -------------------------------------------------------------------------
function cosine_fit

clear all

% -------------------------------------------------------------------------
% Load circadian dataset
%
% Expected format (per *_circ.mat files):
%   column 1 = time-of-day (hours, e.g., 0..24 or 0..23)
%   column 2 = measured hormone concentration
%
% Switch the filename/variable to fit other hormones (P4, LH, FSH).
% -------------------------------------------------------------------------
load E2_circ.mat E2_circ                   % Load the extracted hormonal circadian rhythm data: E2_circ, P4_circ, LH_circ, or FSH_circ   

% Select the hormone time series to fit (here: E2).
% time_data: sampling times in hours; data: hormone values (raw scale).
time_data = E2_circ(:,1);                  % Specify the data: E2_circ, P4_circ, LH_circ, or FSH_circ
data      = E2_circ(:,2);

% Convert to row vectors for convenience.
x1 = (time_data)';    

% Normalize by the mean so the average equals 1 (paper step #2).
% This makes amplitude comparable across hormones and isolates rhythmic pattern.
y1 = (data/mean(data))';

% -------------------------------------------------------------------------
% Define the cosine model to fit (24-hour period)
%
% Model form (hours):
%   y = a1*cos((pi/12)*(x - c1)) + d1
% where:
%   (pi/12) = 2*pi/24 ensures 24-hour periodicity.
%   a1 = amplitude (unitless after normalization)
%   c1 = acrophase (hours; peak timing)
%   d1 = offset (baseline; should be ~1 after normalization)
%
% Paper analogue (days): a*cos(2*pi*(t - b)) + 1.
% -------------------------------------------------------------------------
cosinor_fit = 'a1*cos((pi/12)*(x - c1)) + d1';

% Initial guesses for nonlinear regression [a1, c1, d1].
% c1 initial guess=4 means peak near 4 AM; adjust if fits are unstable.
startpoints1 = [0.5 4 1];

% Fit the cosine curve to normalized data using Curve Fitting Toolbox.
% The resulting parameters can be plugged into the main model's cosine forcing terms.
f1 = fit(x1', y1', cosinor_fit, 'Start', startpoints1)

% Note that the acrophase c1 here  is in hours. 

