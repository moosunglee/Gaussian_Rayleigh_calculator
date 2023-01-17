%% Simple Gaussian optics & Rayleigh approximation
clc, clear;
cd0 = matlab.desktop.editor.getActiveFilename;
cd0 = cd0(1:strfind(cd0,'Run.m')-2);
addpath(genpath(cd0));
%% Every unit: [m, kg, s, W]

% Set default params % Parameters are based on
% https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.103603
    utility = struct;
    % Optics parameters
    utility.lambda = 1064e-9; % [m]
    utility.NA_foc = 0.8;
    utility.NA_det = 0.8;
    utility.Power = 100e-3; 
    % Sample parameters
    utility.radius = 75e-9; % [m]
%     utility.radius = 500e-9; % [m]
    utility.density = 2.196e3;
    utility.RI = GET_RI_SiO2(utility.lambda);
    % Environmental damping parameter
    utility.p = 6.3e2; % [Pa; 1 Bar = 10^5 Bar; 6.3 mBar = 630 Pa]
%     utility.p = 10^-3; % [Pa; 1 Bar = 10^5 Pa; 10^-5 mBar = 10^-3 Pa]
    utility.T = 298; % [K; temperature]

    Calculator = Calc_Gaussian_Rayleigh(utility);
    Calculator.Calc_pGauss('empirical'); % Gaussian optics
    Calculator.Calc_pOptomech;
    % Experimental resonance frequency: 7.54e5, 8.42e5, 2.32e5 
    % (120, 134, 37 kHz in frequency)

    % Backaction frequency: 18.8 kHz * 2 * pi
    Calculator.pOptomech