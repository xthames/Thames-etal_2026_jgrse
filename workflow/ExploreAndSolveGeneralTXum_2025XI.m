
%% An Exploratory Analysis of Earth's Volatile Exchange Rate 
% A scipt that uses a Latin hypercube sampling method to sample Earth's
% mantle's starting parameters to calculate its coupled thermal-water 
% evolution across time.

%% CLOSE and CLEAR and WARNINGS("OFF")
% clear
% close all
% clearvars -except timeDirection viscDep viscStrength modelType varyfR varyphiRum
warning('off')

%% SETTING UP STRUCTS TO HOLD OUR DATA
% struct, ECs for Earth Constants
ECs = [];
% establish the RNG seed
rng('shuffle');
ECs.settings.seed = rng;

% conversion factors
ECs.conversion.degC2K = 273.15;
ECs.conversion.yr2sec = 60*60*24*365.25;
ECs.conversion.sec2yr = 1/ECs.conversion.yr2sec;
ECs.conversion.ppm2massfrac = 1e-6;

%% CONTROLLING WHAT TYPE OF MODEL TO RUN
% time direction for the model set
ECs.settings.timeDirection = timeDirection;
if ~ismember(ECs.settings.timeDirection, ['fwd', 'rev'])
    ECs.settings.timeDirection = 'fwd';
end

% consider both viscosity dependent on water concentration and fugacity
ECs.settings.viscDep = viscDep;
if ~ismember(ECs.settings.viscDep, ['conc', 'fug'])
    ECs.settings.viscDep = 'fug';
end

% strong or weak water-viscosity scaling law
ECs.settings.viscWaterStrength = viscStrength;
if ~ismember(ECs.settings.viscWaterStrength, ['weak', 'strong'])
    ECs.settings.viscWaterStrength = 'strong';
end

% possible variations on the model
ECs.settings.varyfR = 'n';
ECs.settings.varyphiRum = 'n';
ECs.params.W = [1.06 0.14 -0.023 17 2.32];
ECs.settings.varyXp = 'n';

% should we use a baseline model or a variant
ECs.settings.modelType = modelType;
if ~ismember(ECs.settings.modelType, ['base', 'var'])
    % catch for missed input
    ECs.settings.modelType = 'base';
end

% if we're looking at model variants...
if strcmp('var', ECs.settings.modelType)
    % variable regassing fraction/efficiency
    ECs.settings.varyfR = varyfR;
    if ~ismember(ECs.settings.varyfR, ['y', 'n'])
        ECs.settings.varyfR = 'n';
    end

    % fix regassing contribution to upper mantle
    ECs.settings.varyphiRum = varyphiRum;
    if ~ismember(ECs.settings.varyphiRum, ['y', 'n'])
        ECs.settings.varyphiRum = 'n';
    end
end

% display what type of model will be run
SectionBreak()
disp(['SETTINGS for GENERAL ', upper(ECs.settings.modelType),' ',upper(ECs.settings.viscDep),' models:'])
disp(['-- Water-Viscosity Strength: ', upper(ECs.settings.viscWaterStrength)])
disp(['-- Variable fR: ', upper(ECs.settings.varyfR)])
disp(['-- Variable phiRum: ', upper(ECs.settings.varyphiRum)])

%% NUMBER OF MODELS TO RUN
ECs.MC.iterations = 1e5;
SectionBreak()
disp(['RUNNING: ',num2str(ECs.MC.iterations),' scenarios...'])

%% ESTIMATING HOW LONG THE CODE WILL TAKE TO RUN...
% time estimator
ECs.ETC.time = 3.5e-2 * exp(2.207*log10(ECs.MC.iterations));

% conditional
if ECs.ETC.time < 60
    ECs.ETC.value = ECs.ETC.time;
    ECs.ETC.unit = 'seconds';
elseif ECs.ETC.time >= 60 && ECs.ETC.time < 60*60
    ECs.ETC.value = ECs.ETC.time / 60;
    ECs.ETC.unit = 'minutes';
elseif ECs.ETC.time >= 60*60 && ECs.ETC.time < 60*60*24
    ECs.ETC.value = ECs.ETC.time / (60*60);
    ECs.ETC.unit = 'hours';
else
    ECs.ETC.value = ECs.ETC.time / (60*60*24);
    ECs.ETC.unit = 'days';
end

% display
disp(['ESTIMATED TIME: ', num2str(ECs.ETC.value), ' ', ECs.ETC.unit])

%% AVERAGE EARTH CONSTANTS
ECs.avg.rho = 4000;                                                                     % average bulk mantle density [kg m^-3]
ECs.avg.k = 3.2;                                                                        % thermal conductivity [W m^-1 K^-1]
ECs.avg.c_p = 1200;                                                                     % specific heat [J kg^-1 K^-1]
ECs.avg.kappa = 7.6e-7;                                                                 % thermal diffusivity [m^2 s^-1]
ECs.avg.alpha = 2.5e-5;                                                                 % thermal expension [K^-1]
ECs.avg.g = 9.81;                                                                       % gravitational acceleration [m s^-2]
ECs.avg.Rg = 8.31;                                                                      % Gas Constant [J mol^-1 K^-1]
ECs.avg.M_E = 5.9722e24;                                                                % Earth's mass [kg]
ECs.avg.R_E = 6371e3;                                                                   % Earth's radius [m]
ECs.avg.R_Ec = 3491e3;                                                                  % Earth's core radius [m]
ECs.avg.R_Em = ECs.avg.R_E - ECs.avg.R_Ec;                                              % Earth's mantle thickess [m]
ECs.avg.h_um = 410e3;                                                                   % depth to bottom of upper mantle [m]
ECs.avg.h_tz = 660e3;                                                                   % depth to bottom of transition zone [m]
ECs.avg.A_O = 0.71*4*pi*ECs.avg.R_E^2;                                                  % Earth's oceanic basin surface area [m^2]
ECs.avg.V_M = 4/3*pi*(ECs.avg.R_E^3 - ECs.avg.R_Ec^3);                                  % Earth's mantle volume [m^3]
ECs.avg.Vum = 4/3*pi*(ECs.avg.R_E^3 - (ECs.avg.R_E - ECs.avg.h_um)^3);                  % volume of upper mantle [m]
ECs.avg.Vtz = 4/3*pi*((ECs.avg.R_E - ECs.avg.h_um)^3 - (ECs.avg.R_E - ECs.avg.h_tz)^3); % volume of transition zone [m] 
ECs.avg.Vlm = 4/3*pi*((ECs.avg.R_E - ECs.avg.h_tz)^3 - ECs.avg.R_Ec^3);                 % volume of lower mantle [m]
ECs.avg.L = 6000e3;                                                                     % average plate width [m] 
ECs.avg.T = mean([38500e3 56800e3]);                                                    % approx. length of subduction trenches [m, van Keken (2011) says 38,500km, Lallemand (1999) says 56,800km]
ECs.avg.rhop = 3315;                                                                    % density of subducting plate, [kg m^-3] (Alfonso et al. 2007)
ECs.avg.rhoum = 3300;									% density of upper mantle at approx surface pressure [kg m^-3]

%% VISCOSITY PARAMETERS
ECs.visc.Ra_c = 1e3;                            % the critical Rayleigh number, below which convection does not happen in Earth's mantle [-]
ECs.visc.Xc = 1*ECs.conversion.ppm2massfrac;    % normalization conc. [kg/kg], Fei et al. (2013)
ECs.visc.rW = 1/3;                              % weak water visc exponent, eta(T,Xum) [-]
ECs.visc.rS = 1;                                % strong water visc exponent, eta(T,Xum) [-]
if strcmp('strong', ECs.settings.viscWaterStrength)
    ECs.visc.r = ECs.visc.rS;
else
    ECs.visc.r = ECs.visc.rW;
end
ECs.visc.fugC0 = -8.0;                          % water fugacity term, zeroth order [-]
ECs.visc.fugC1 = 4.4;                           % water fugacity term, first order [-]
ECs.visc.fugC2 = -0.57;                         % water fugacity term, second order [-]
ECs.visc.fugC3 = 0.033;                         % water fugacity term, third order [-]
ECs.visc.ppm2COH = 16.3;                        % convert between wt ppm H2O to H / 10^6 Si [-]
ECs.params.beta = 1/3;                          % the beta parameter, defining the relationship between Nusselt and Rayleigh numbers

%% CONSERVATION OF ENERGY PARAMETERS (from Treatise on Geophysics)
% concentration [ppm]
ECs.consE.C_Us = [.020 .018 .021 .020 .022 .017 .014]; ECs.consE.U235pct = 0.72e-2; ECs.consE.U238pct = 99.27e-2;
ECs.consE.C_Ths = [.069 .060 .079 .079 .083 .063 .055]; ECs.consE.Th232pct = 100e-2;
ECs.consE.C_Ks = [270 217 264 240 261 190 166]; ECs.consE.K40pct = 0.0117e-2;

% heating rate [W kg^-1]
ECs.consE.H_U238 = 95.13e-6;
ECs.consE.H_U235 = 568.47e-6;
ECs.consE.H_Th232 = 26.3e-6;
ECs.consE.H_K40 = 28.47e-6;

% half-lives [s]
ECs.consE.tau_U238 = 4.46e9 * ECs.conversion.yr2sec;
ECs.consE.tau_U235 = 0.704e9 * ECs.conversion.yr2sec;
ECs.consE.tau_Th232 = 14.0e9 * ECs.conversion.yr2sec;
ECs.consE.tau_K40 = 1.26e9 * ECs.conversion.yr2sec;

% ranges for fraction of heat production in bulk mantle and crust [TW]
ECs.consE.heatbudget_convectingmantle = [9 17];
ECs.consE.heatbudget_crust = [7 8];

%% DEGASSING PARAMETERS
% degassing efficiency [-]
ECs.D.fD = 1;

% melt depth constants
ECs.D.z1 = 286;                                     % [m (deg C)^-1]
ECs.D.z2 = 164/ECs.conversion.ppm2massfrac;         % [m (kg/kg)^-1]
ECs.D.z3 = -3.2e5;                                  % [m]

%% PRESENT-DAY GENERAL
% time
ECs.pd.ageEarth = 4.543e9*ECs.conversion.yr2sec;                    % [s]
ECs.pd.t = 4.543e9*ECs.conversion.yr2sec;                           % [s]

% temperature and surface ocean mass
ECs.pd.T = 1350 + ECs.conversion.degC2K;                            % [K]
ECs.pd.Tsurf = 15 + ECs.conversion.degC2K;                          % [K]
ECs.pd.oceanmass = 1.39e21;                                         % [kg]

% viscosity, plate speed, heat flux
ECs.pd.eta_um = 1e21;                                               % [Pa s]
ECs.pd.u = 4.5 * ECs.conversion.sec2yr * 0.01;                      % [m s^-1]
ECs.pd.Qs = 32e12;                                                  % [W], note that 32TW from Jaupart (2015) uses an area of 300 million km^2

% scaling constants for eta, u, Qs based on present-day values
ECs.pd.delta = 2 * sqrt((ECs.avg.kappa * ECs.avg.L) / ECs.pd.u);    % [m]
ECs.pd.E = 300e3;                                                   % [J mol^-3]
ECs.pd.Ra = (ECs.avg.alpha * ECs.avg.rho * ECs.avg.g * (ECs.pd.T - ECs.pd.Tsurf) * (ECs.avg.R_Em^3)) / (ECs.avg.kappa * ECs.pd.eta_um);
ECs.pd.uFromReference = (ECs.avg.kappa / ECs.avg.R_Em) * ((ECs.avg.L / (pi * ECs.avg.R_Em))^(ECs.params.beta)) * ECs.pd.Ra^(2 * ECs.params.beta);
ECs.scale.u0 = ECs.pd.u / ECs.pd.uFromReference;
ECs.pd.QsFromReference = 2 * ECs.avg.A_O * ECs.avg.k * (ECs.pd.T - ECs.pd.Tsurf) * sqrt(ECs.pd.u / (pi * ECs.avg.L * ECs.avg.kappa));
ECs.scale.Qs0 = ECs.pd.Qs / ECs.pd.QsFromReference;

%% PRESENT-DAY WATER DISTRIBUTION IN THE MANTLE
% water distributions in the mantle layers
% -- UM: ~100 ppm, can be 1-64 ppm but as high as 300ppm. ~0.04 OM
% -- TZ: using Ohtani (2019). somewhere between 0.2 and 1 OM
% -- LM: <<0.1 wt%, Hirschmann (2006) says 20 ppm. <2 OM. 
ECs.pd.H2Oconcs_tz = [0.1 0.3 0.01 0.1 0.1 0.2 2 0.5 1 0.1 1].*1e-2;    % [wt%], from Ohtani (2019)
ECs.pd.Xtz = median(ECs.pd.H2Oconcs_tz);
ECs.pd.H2Oconcs_lm = [20*ECs.conversion.ppm2massfrac 0.1e-2];        % [kg/kg], Hirschmann (2006), [upper bound] from Ohtani (2019)
ECs.pd.Xlm = median(ECs.pd.H2Oconcs_lm);
ECs.pd.Xum = 100 * ECs.conversion.ppm2massfrac;

% water concentration in mantle layers, in [kg/kg] or [kg]
ECs.pd.MH2O_um = ECs.pd.Xum * ECs.avg.rho * ECs.avg.Vum;
ECs.pd.MH2O_tz = ECs.pd.Xtz * ECs.avg.rho * ECs.avg.Vtz;
ECs.pd.MH2O_lm = ECs.pd.Xlm * ECs.avg.rho * ECs.avg.Vlm;
ECs.pd.MH2O_totalM = ECs.pd.MH2O_um + ECs.pd.MH2O_tz + ECs.pd.MH2O_lm;

% fraction of water stored in upper mantle, present-day concentration and fugacity
ECs.pd.MH2Ofrac_um = ECs.pd.MH2O_um / ECs.pd.MH2O_totalM;
ECs.pd.H2Oconc = ECs.pd.Xum / ECs.visc.Xc;
ECs.scale.eta0_conc = ECs.pd.eta_um / (((ECs.pd.H2Oconc)^(-ECs.visc.r))*exp(ECs.pd.E/(ECs.avg.Rg*ECs.pd.T)));
ECs.pd.H2Ofug = exp(ECs.visc.fugC0 + ...
                    ECs.visc.fugC1 * log(ECs.visc.ppm2COH * ECs.pd.Xum / ECs.conversion.ppm2massfrac) + ...
                    ECs.visc.fugC2 * (log(ECs.visc.ppm2COH * ECs.pd.Xum / ECs.conversion.ppm2massfrac))^2 + ...
                    ECs.visc.fugC3 * (log(ECs.visc.ppm2COH * ECs.pd.Xum / ECs.conversion.ppm2massfrac))^3);
ECs.scale.eta0_fug = ECs.pd.eta_um / (((ECs.pd.H2Ofug)^(-ECs.visc.r))*exp(ECs.pd.E/(ECs.avg.Rg*ECs.pd.T)));

%% RANGES FOR THE LATIN HYPERCUBE
% time [s]
if strcmp('rev', ECs.settings.timeDirection)
    ECs.range.t = [ECs.pd.t 0];
else
    ECs.range.t = [0 ECs.pd.t];
end

% temperature [K]
ECs.range.Tpd = [ECs.pd.T-80 ECs.pd.T+95];
ECs.range.T0 = [min(ECs.range.Tpd - ECs.conversion.degC2K) 1650] + ECs.conversion.degC2K;

% water
ECs.range.totalH2Omass = [1.25 10] .* ECs.pd.oceanmass;
ECs.range.H2Opartition = [0.001 0.999];
ECs.range.Xp = [860 7800] .* ECs.conversion.ppm2massfrac; % when reduced to hydration depth, could be [0.2 - 10]%
ECs.range.fum = [0.35 63.3] ./ 100;

% variable regassing [-], bounds from Magni (2014)
ECs.range.fR = [0.001 0.999];
if ECs.settings.varyfR == 'y'
    ECs.range.fR0 = [1 3]; 
    ECs.range.MagniCLinearScale = [0.6 1];
end

% variable regassing contribution to upper mantle
ECs.range.phiRum = [0.001 0.999];

% radiogenics [kg/kg]
ECs.range.fRGinMantle = [(0.9*min(ECs.consE.heatbudget_convectingmantle))/(0.9*min(ECs.consE.heatbudget_convectingmantle) + 1.1*max(ECs.consE.heatbudget_crust)), ...
                         (1.1*max(ECs.consE.heatbudget_convectingmantle))/(1.1*max(ECs.consE.heatbudget_convectingmantle) + 0.9*min(ECs.consE.heatbudget_crust))];
ECs.range.C_U235 = [min(ECs.consE.C_Us)-0.1*min(ECs.consE.C_Us) max(ECs.consE.C_Us)+0.1*max(ECs.consE.C_Us)] .* ECs.consE.U235pct .* ECs.conversion.ppm2massfrac;
ECs.range.C_U238 = [min(ECs.consE.C_Us)-0.1*min(ECs.consE.C_Us) max(ECs.consE.C_Us)+0.1*max(ECs.consE.C_Us)] .* ECs.consE.U238pct .* ECs.conversion.ppm2massfrac;
ECs.range.C_Th232 = [min(ECs.consE.C_Ths)-0.1*min(ECs.consE.C_Ths) max(ECs.consE.C_Ths)+0.1*max(ECs.consE.C_Ths)] .* ECs.consE.Th232pct .* ECs.conversion.ppm2massfrac;
ECs.range.C_K40 = [min(ECs.consE.C_Ks)-0.1*min(ECs.consE.C_Ks) max(ECs.consE.C_Ks)+0.1*max(ECs.consE.C_Ks)] .* ECs.consE.K40pct .* ECs.conversion.ppm2massfrac;

% activation energy [J mol^-3]
ECs.range.E = [260 450] .* 1e3;

% present-day upper mantle viscosity (10^[this value], Pa s)
ECs.range.eta_pd = [20 22];

%% LATIN HYPERCUBE SAMPLING FOR PARAMETERS
disp('ESTABLISHING LHS...')
% design the LHS
if ECs.settings.varyfR == 'y'
    ECs.lhs.design = lhsdesign(ECs.MC.iterations, 14);
else
    ECs.lhs.design = lhsdesign(ECs.MC.iterations, 13);
end

% fill with parameters
ECs.lhs.params = NaN(size(ECs.lhs.design));
for i = 1:size(ECs.lhs.design, 1)
    % temperature
    if strcmp('rev', ECs.settings.timeDirection)
        ECs.lhs.params(i, 1) = interp1([0, 1], ECs.range.Tpd, ECs.lhs.design(i, 1));
    else
        ECs.lhs.params(i, 1) = interp1([0, 1], ECs.range.T0, ECs.lhs.design(i, 1));
    end
    % total H2O mass
    ECs.lhs.params(i, 2) = interp1([0, 1], ECs.range.totalH2Omass, ECs.lhs.design(i, 2));
    % H2O partition between mantle and surface
    ECs.lhs.params(i, 3) = interp1([0, 1], ECs.range.H2Opartition, ECs.lhs.design(i, 3));
    % water concentration on the plate Xp
    ECs.lhs.params(i, 4) = interp1([0, 1], ECs.range.Xp, ECs.lhs.design(i, 4));
    % fraction of radiogenics in the mantle
    ECs.lhs.params(i, 5) = interp1([0, 1], ECs.range.fRGinMantle, ECs.lhs.design(i, 5));
    % radiogenic concentrations
    ECs.lhs.params(i, 6) = interp1([0, 1], ECs.range.C_U235, ECs.lhs.design(i, 6));
    ECs.lhs.params(i, 7) = interp1([0, 1], ECs.range.C_U238, ECs.lhs.design(i, 7));
    ECs.lhs.params(i, 8) = interp1([0, 1], ECs.range.C_Th232, ECs.lhs.design(i, 8));
    ECs.lhs.params(i, 9) = interp1([0, 1], ECs.range.C_K40, ECs.lhs.design(i, 9));
    % activation energy
    ECs.lhs.params(i, 10) = interp1([0, 1], ECs.range.E, ECs.lhs.design(i, 10));
    % present-day viscosity
    ECs.lhs.params(i, 11) = interp1([0, 1], ECs.range.eta_pd, ECs.lhs.design(i, 11));
    % regassing efficiency, frac (contribution of R to) upper mantle
    if ECs.settings.varyfR == 'y'
	ECs.lhs.params(i, 12) = interp1([0, 1], ECs.range.fR0, ECs.lhs.design(i, 12));
        ECs.lhs.params(i, 13) = interp1([0, 1], ECs.range.MagniCLinearScale, ECs.lhs.design(i, 13));
	if ECs.settings.varyphiRum == 'y'
            ECs.lhs.params(i, 14) = interp1([0, 1], ECs.range.phiRum, ECs.lhs.design(i, 14));
	else
            ECs.lhs.params(i, 14) = interp1([0, 1], ECs.range.fum, ECs.lhs.design(i, 14));
	end
    else
        ECs.lhs.params(i, 12) = interp1([0, 1], ECs.range.fR, ECs.lhs.design(i, 12));
        ECs.lhs.params(i, 13) = interp1([0, 1], ECs.range.fum, ECs.lhs.design(i, 13));
    end
end

%% SOLVE THE CONSERVATION EQUATIONS
disp('SOLVING CONSERVATION EQUATIONS...')
% output structs for the parfor loop
solved(ECs.MC.iterations) = struct();

tic
parfor i = 1:ECs.MC.iterations
    % set this loop's constants to our global set, clear out solved time series
    LCs = ECs;
    LTS = [];
    
    % --- temperature
    LCs.params.T = LCs.lhs.params(i, 1);

    % --- water
    LCs.params.totalH2Omass = LCs.lhs.params(i, 2);
    
    % --- water concentration in the plate OR hydrated penetration depth
    LCs.params.Xp = LCs.lhs.params(i, 4);
    
    % --- fraction of RG in mantle
    LCs.params.fRGinMantle = LCs.lhs.params(i, 5);
    
    % --- radiogenics
    LCs.params.C_U235 = LCs.params.fRGinMantle * LCs.lhs.params(i, 6);
    LCs.params.C_U238 = LCs.params.fRGinMantle * LCs.lhs.params(i, 7);
    LCs.params.C_Th232 = LCs.params.fRGinMantle * LCs.lhs.params(i, 8);
    LCs.params.C_K40 = LCs.params.fRGinMantle * LCs.lhs.params(i, 9);
    
    % --- activation energy
    LCs.params.E = LCs.lhs.params(i, 10);
    
    % --- upper mantle viscosity
    LCs.params.eta_pd = 10^(LCs.lhs.params(i, 11));
    
    % --- volatile
    if LCs.settings.varyfR == 'y'
        LCs.scale.fR0 = 10^(LCs.lhs.params(i, 12));
        LCs.scale.MagniCScale = LCs.lhs.params(i, 13);
	if LCs.settings.varyphiRum == 'y'
	    LCs.params.phiRum = LCs.lhs.params(i, 14);
        else
	    LCs.params.fum = LCs.lhs.params(i, 14);
	end
    else
        LCs.params.fR = LCs.lhs.params(i, 12);
        LCs.params.fum = LCs.lhs.params(i, 13);
    end
    
    % -- time dependent
    if strcmp('fwd', LCs.settings.timeDirection)
        % water
        LCs.params.H2Opartition = LCs.lhs.params(i, 3);
        if LCs.settings.varyphiRum == 'y'
            LCs.params.Xum = LCs.params.totalH2Omass * LCs.params.H2Opartition * (LCs.avg.Vum/LCs.avg.V_M) / (LCs.avg.rhoum * LCs.avg.Vum);
        else
            LCs.params.Xum = LCs.params.totalH2Omass * LCs.params.H2Opartition * LCs.params.fum / (LCs.avg.rhoum * LCs.avg.Vum);
        end

        % scaling
        if strcmp('conc', LCs.settings.viscDep)
            LCs.scale.eta0 = LCs.params.eta_pd / (((LCs.pd.H2Oconc)^(-LCs.visc.r))*exp(LCs.params.E / (LCs.avg.Rg * LCs.pd.T)));
        else
            LCs.scale.eta0 = LCs.params.eta_pd / (((LCs.pd.H2Ofug)^(-LCs.visc.r))*exp(LCs.params.E / (LCs.avg.Rg * LCs.pd.T)));
        end
        LCs.scale.Ra = (LCs.avg.alpha * LCs.avg.rho * LCs.avg.g * (LCs.pd.T - LCs.pd.Tsurf) * LCs.avg.R_Em^3) / (LCs.avg.kappa * LCs.params.eta_pd);
        LCs.scale.uFromReference = (LCs.avg.kappa / LCs.avg.R_Em) * ((LCs.avg.L / (pi * LCs.avg.R_Em))^(LCs.params.beta)) * LCs.scale.Ra^(2 * LCs.params.beta);
        LCs.scale.u = LCs.pd.u / LCs.scale.uFromReference;
        LCs.scale.Qs = LCs.scale.Qs0;
    else
        % water
        LCs.params.H2Opartition = (LCs.params.totalH2Omass - LCs.pd.oceanmass) / LCs.params.totalH2Omass;
        LCs.params.Xum = (LCs.params.totalH2Omass - LCs.pd.oceanmass) * LCs.params.fum / (LCs.avg.rhoum * LCs.avg.Vum);

        % scaling
        if strcmp('conc', LCs.settings.viscDep)
            LCs.scale.H2Oconc = LCs.params.Xum / LCs.visc.Xc;
            LCs.scale.eta0 = LCs.params.eta_pd / (((LCs.scale.H2Oconc)^(-LCs.visc.r))*exp(LCs.params.E/(LCs.avg.Rg*LCs.params.T)));
        else
            LCs.scale.H2Ofug = exp(LCs.visc.fugC0 + ...
                                   LCs.visc.fugC1 * log(LCs.visc.ppm2COH * LCs.params.Xum / LCs.conversion.ppm2massfrac) + ...
                                   LCs.visc.fugC2 * (log(LCs.visc.ppm2COH * LCs.params.Xum / LCs.conversion.ppm2massfrac))^2 + ...
                                   LCs.visc.fugC3 * (log(LCs.visc.ppm2COH * LCs.params.Xum / LCs.conversion.ppm2massfrac))^3);
            LCs.scale.eta0 = LCs.params.eta_pd / (((LCs.scale.H2Ofug)^(-LCs.visc.r))*exp(LCs.params.E/(LCs.avg.Rg*LCs.params.T)));
        end
        LCs.scale.Ra = (LCs.avg.alpha * LCs.avg.rho * LCs.avg.g * (LCs.params.T - LCs.pd.Tsurf) * LCs.avg.R_Em^3) / (LCs.avg.kappa * LCs.params.eta_pd);
        LCs.scale.uFromReference = (LCs.avg.kappa / LCs.avg.R_Em) * ((LCs.avg.L / (pi * LCs.avg.R_Em))^(LCs.params.beta)) * LCs.scale.Ra^(2 * LCs.params.beta);
        LCs.scale.u = LCs.pd.u / LCs.scale.uFromReference;
        LCs.scale.QsFromReference = 2 * LCs.avg.A_O * LCs.avg.k * (LCs.params.T - LCs.pd.Tsurf) * sqrt(LCs.pd.u/(pi*LCs.avg.L*LCs.avg.kappa));
        LCs.scale.Qs = LCs.pd.Qs / LCs.scale.QsFromReference;
    end
    
    % initial conditions, ode solver options
    LCs.ICs = [LCs.params.T, LCs.params.Xum];
    LCs.odeOpts = odeset('RelTol', 1e-10, 'NonNegative', [1 2]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % solve the complete time series for T and Xum
    [LCs.TS.t, Y] = ode15s(@(t, Y) ConservationEquations(t, Y, LCs), LCs.range.t, LCs.ICs, LCs.odeOpts);
    LCs.TS.T = Y(:, 1); LCs.TS.Xum = Y(:, 2);

    % MELT DEPTH
    LCs.TS.z = LCs.D.z1.*LCs.TS.T + LCs.D.z2.*LCs.TS.Xum + LCs.D.z3;

    % HEAT PRODUCTION
    LCs.TS.H = LCs.avg.rho .* LCs.avg.V_M .* ...
              (LCs.params.C_U238 .* LCs.consE.H_U238 .* exp(log(2) .* (LCs.pd.ageEarth - LCs.TS.t) ./ LCs.consE.tau_U238) + ...
               LCs.params.C_U235 .* LCs.consE.H_U235 .* exp(log(2) .* (LCs.pd.ageEarth - LCs.TS.t) ./ LCs.consE.tau_U235) + ...
               LCs.params.C_Th232 .* LCs.consE.H_Th232 .* exp(log(2) .* (LCs.pd.ageEarth - LCs.TS.t) ./ LCs.consE.tau_Th232) + ...
               LCs.params.C_K40 .* LCs.consE.H_K40 .* exp(log(2) .* (LCs.pd.ageEarth - LCs.TS.t) ./ LCs.consE.tau_K40));
    
    % VISCOSITY
    if strcmp('conc', LCs.settings.viscDep)
        LCs.TS.H2Oconc = LCs.TS.Xum ./ LCs.visc.Xc;
        LCs.TS.eta_um = LCs.scale.eta0 .* (((LCs.TS.H2Oconc).^(-LCs.visc.r)) .* exp(LCs.params.E ./ (LCs.avg.Rg .* LCs.TS.T)));
    else
        LCs.TS.H2Ofug = exp(LCs.visc.fugC0 + ...
                            LCs.visc.fugC1 .* log(LCs.visc.ppm2COH .* LCs.TS.Xum ./ LCs.conversion.ppm2massfrac) + ...
                            LCs.visc.fugC2 .* (log(LCs.visc.ppm2COH .* LCs.TS.Xum ./ LCs.conversion.ppm2massfrac)).^2 + ...
                            LCs.visc.fugC3 .* (log(LCs.visc.ppm2COH .* LCs.TS.Xum ./ LCs.conversion.ppm2massfrac)).^3);
        LCs.TS.eta_um = LCs.scale.eta0 .* (((LCs.TS.H2Ofug).^(-LCs.visc.r)) .* exp(LCs.params.E ./ (LCs.avg.Rg .* LCs.TS.T)));
    end
    
    % RAYLEIGH NUMBER
    LCs.TS.Ra = ((LCs.avg.alpha .* LCs.avg.rho .* LCs.avg.g) .* (LCs.TS.T - LCs.pd.Tsurf) .* LCs.avg.R_Em.^3) ./ (LCs.avg.kappa .* LCs.TS.eta_um);

    % PLATE SPEED
    LCs.TS.u = LCs.scale.u .* (LCs.avg.kappa ./ LCs.avg.R_Em) .* (LCs.avg.L ./ (pi .* LCs.avg.R_Em)).^(LCs.params.beta) .* LCs.TS.Ra.^(2.*LCs.params.beta);

    % PLATE THICKNESS
    LCs.TS.delta = 2 .* sqrt((LCs.avg.kappa .* LCs.avg.L) ./ LCs.TS.u);

    % TOTAL HEAT FLOW THROUGH SURFACE AREA
    LCs.TS.Qs = 2 .* LCs.scale.Qs .* LCs.avg.A_O .* LCs.avg.k .* (LCs.TS.T - LCs.pd.Tsurf) .* sqrt(LCs.TS.u ./ (pi .* LCs.avg.L .* LCs.avg.kappa));

    % DEGASSING EXCHANGE RATE
    LCs.TS.D = LCs.avg.A_O .* LCs.D.fD .* (LCs.TS.z ./ LCs.avg.L) .* LCs.TS.u .* LCs.avg.rhoum .* LCs.TS.Xum;
 
    % PLATE WATER CONCENTRATION
    LCs.TS.Xp = LCs.params.Xp .* ones(size(LCs.TS.t));
    
    % REGASSING EXCHANGE RATE
    if LCs.settings.varyfR == 'y'
    	% -- variant: variable fR
        a = ((LCs.TS.delta ./ LCs.params.W(5)).^2) ./ LCs.avg.kappa .* LCs.conversion.sec2yr .* 1e-6;
        W = 1e5 .* ((LCs.params.W(1).*(LCs.TS.u .* LCs.conversion.yr2sec .* 100)) + ...
                    (LCs.params.W(2).*a) + ...
                    (LCs.params.W(3).*(LCs.TS.T - LCs.conversion.degC2K)) + ...
                    LCs.params.W(4).*LCs.scale.MagniCScale);
        W0 = LCs.TS.Xp .* LCs.avg.rhop .* (LCs.avg.A_O ./ LCs.avg.T);
        LCs.TS.fR = LCs.scale.fR0 * (W ./ W0);
        LCs.TS.fR(LCs.TS.fR <= 0) = eps; LCs.TS.fR(LCs.TS.fR >= 1) = 1-eps;
    else
        LCs.TS.fR = LCs.params.fR .* ones(size(LCs.TS.t));
    end
    LCs.TS.R = LCs.avg.A_O .* LCs.TS.fR .* (LCs.TS.delta ./ LCs.avg.L) .* LCs.TS.u .* LCs.avg.rhop .* LCs.TS.Xp;
    
    % OCEAN MASS
    if LCs.settings.varyphiRum == 'y'
    	% -- variant: variable phiRum
        LCs.TS.phiRum = LCs.params.phiRum .* ones(size(LCs.TS.R));
        deltat = gradient(LCs.TS.t);
        LCs.TS.mantleH2Omass = NaN(size(LCs.TS.t));
        for k = 1:length(LCs.TS.mantleH2Omass)
            if k == 1
                LCs.TS.mantleH2Omass(k) = (LCs.params.totalH2Omass * LCs.params.H2Opartition) + (LCs.TS.R(k) - LCs.TS.D(k)).*deltat(k);
            else
                LCs.TS.mantleH2Omass(k) = LCs.TS.mantleH2Omass(k-1) + (LCs.TS.R(k) - LCs.TS.D(k)).*deltat(k);
            end
        end
    else
        LCs.TS.mantleH2Omass = (LCs.TS.Xum .* LCs.avg.rhoum .* LCs.avg.Vum ./ LCs.params.fum);
    end
    LCs.TS.oceanmass = LCs.params.totalH2Omass - LCs.TS.mantleH2Omass;

    % UREY RATIO
    LCs.TS.UreyRatio = LCs.TS.H ./ LCs.TS.Qs; 
    
    % FRACTION OF TOTAL MANTLE WATER MASS IN UPPER MANTLE
    LCs.TS.fH2Oum = (LCs.TS.Xum .* LCs.avg.rhoum .* LCs.avg.Vum) ./ (LCs.params.totalH2Omass - LCs.TS.oceanmass);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set all the timeseries data into solved(i).TS, with human units
    % -- time check, did it run the full simulation?
    if (LCs.TS.t(1) == LCs.range.t(1)) && (LCs.TS.t(end) == LCs.range.t(end))
        % -- make sure Xum, oceanmass, Ra are always positive and above
        % -- thresholds, as well as that always (mantle + surface = total water)
        if all(LCs.TS.Xum > 0) && ...
           all(LCs.TS.oceanmass > 0) && all(ismembertol(LCs.TS.oceanmass + LCs.TS.mantleH2Omass, LCs.params.totalH2Omass)) && ...
           all(LCs.TS.Ra - LCs.visc.Ra_c > 0)
            solved(i).MC.fullSim = true;
        else
            solved(i).MC.fullSim = false;
        end
    else
        solved(i).MC.fullSim = false;
    end

    % -- time [Gyr]
    solved(i).TS.t = LCs.TS.t .* LCs.conversion.sec2yr ./ 1e9;

    % -- temperature [degC]
    solved(i).TS.T = LCs.TS.T - LCs.conversion.degC2K;
    
    % -- upper mantle water concentration
    solved(i).TS.Xum = LCs.TS.Xum ./ LCs.conversion.ppm2massfrac;

    % -- melt depth [km]
    solved(i).TS.z = LCs.TS.z ./ 1e3;

    % -- total heat production [TW]
    solved(i).TS.H = LCs.TS.H ./ 1e12;

    % -- viscosity [Pa s]
    solved(i).TS.eta_um = LCs.TS.eta_um;

    % -- Rayleigh Number [-]
    solved(i).TS.Ra = LCs.TS.Ra;

    % -- plate speed [cm/yr]
    solved(i).TS.u = LCs.TS.u ./ LCs.conversion.sec2yr .* 100;

    % -- total heat flow through surface [TW]
    solved(i).TS.Qs = LCs.TS.Qs ./ 1e12;

    % volatile movement [kg/yr]
    solved(i).TS.D = LCs.TS.D ./ LCs.conversion.sec2yr;
    solved(i).TS.R = LCs.TS.R ./ LCs.conversion.sec2yr;
    solved(i).TS.RminusD = (solved(i).TS.R - solved(i).TS.D);
    
    % Urey Ratio [-]
    solved(i).TS.UreyRatio = LCs.TS.UreyRatio;
    
    % model variant: penetration depth
    % -- plate thickness [km]
    solved(i).TS.delta = LCs.TS.delta ./ 1e3;

    % water masses [OM]
    solved(i).TS.mantleH2Omass = LCs.TS.mantleH2Omass ./ LCs.pd.oceanmass;
    solved(i).TS.oceanmass = LCs.TS.oceanmass ./ LCs.pd.oceanmass;
    
    % model variant: variable Xp
    % -- plate thickness [km]
    solved(i).TS.Xp = LCs.TS.Xp ./ LCs.conversion.ppm2massfrac;
    
    % model variant: variable fR
    % -- variable fR [%]
    solved(i).TS.fR = LCs.TS.fR * 100;

    % model variant: variable phiRum
    % -- variable phiRum [%]
    solved(i).TS.fH2Oum = LCs.TS.fH2Oum .* 100;
    if LCs.settings.varyphiRum == 'y'
        solved(i).TS.phiRum = LCs.TS.phiRum * 100;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set all the params into solved(i).MC, with human units
    % -- total ocean mass [OM]
    solved(i).MC.totalH2Omass = LCs.params.totalH2Omass / LCs.pd.oceanmass;

    % -- partition of total H2O mass in mantle [%]
    solved(i).MC.H2Opartition = LCs.params.H2Opartition * 100;

    % -- fraction of RG in the mantle of total [%]
    solved(i).MC.fRGinMantle = LCs.params.fRGinMantle * 100;
    
    % -- concentration of RGs at present-day [ppb]
    solved(i).MC.C_U235 = LCs.params.C_U235 / LCs.conversion.ppm2massfrac * 1e3;
    solved(i).MC.C_U238 = LCs.params.C_U238 / LCs.conversion.ppm2massfrac * 1e3;
    solved(i).MC.C_Th232 = LCs.params.C_Th232 / LCs.conversion.ppm2massfrac * 1e3;
    solved(i).MC.C_K40 = LCs.params.C_K40 / LCs.conversion.ppm2massfrac * 1e3;
    
    % -- (present-day) activation energy [kJ mol^3]
    solved(i).MC.E = LCs.params.E / 1e3;

    % -- selected present-day viscosity [Pa s]
    solved(i).MC.eta_pd = LCs.params.eta_pd;
    
    % -- scaling parameters
    solved(i).MC.eta0 = LCs.scale.eta0;
    solved(i).MC.u_scale = LCs.scale.u;
    solved(i).MC.Qs_scale = LCs.scale.Qs;
    if LCs.settings.varyfR == 'y'
        solved(i).MC.fR0 = LCs.scale.fR0;
        solved(i).MC.MagniCScale = LCs.scale.MagniCScale;
	if LCs.settings.varyphiRum == 'y'
	    solved(i).MC.phiRum = LCs.params.phiRum * 100;
        else
	    solved(i).MC.fum = LCs.params.fum * 100;
	end
    else
        solved(i).MC.fR = LCs.params.fR * 100;
        solved(i).MC.fum = LCs.params.fum * 100;
    end

    % -- time when R - D = 0, geochemistry says must transition from D to R in forward direction
    if strcmp('fwd', LCs.settings.timeDirection) && solved(i).MC.fullSim && any(solved(i).TS.RminusD < 0) && any(solved(i).TS.RminusD > 0)
        % spline fit, find the approximate transition time
        LCs.idx.RminusDfit = fit(solved(i).TS.t, solved(i).TS.RminusD, 'smoothingspline');
        if max(LCs.idx.RminusDfit(solved(i).TS.t)) < 0 || min(LCs.idx.RminusDfit(solved(i).TS.t)) > 0
            solved(i).MC.tRminusD0 = NaN;
        end
        if ~isfield(solved(i).MC, 'tRminusD0')
            % index of largest negative, next smallest positive values for volatile exchange
            LCs.idx.mostnegRminusDval = min(solved(i).TS.RminusD(solved(i).TS.RminusD < 0));
            LCs.idx.mostnegRminusD = find(solved(i).TS.RminusD == LCs.idx.mostnegRminusDval);
            LCs.idx.neg2posRminusD = solved(i).TS.RminusD(LCs.idx.mostnegRminusD:end);
            if ~any(LCs.idx.neg2posRminusD > 0)
                solved(i).MC.tRminusD0 = fzero(@(t) LCs.idx.RminusDfit(t), [solved(i).TS.t(1), solved(i).TS.t(LCs.idx.mostnegRminusD-1)], optimset('Display','off'));
            end
            if ~isfield(solved(i).MC, 'tRminusD0')
                LCs.idx.nextleastposRminusDval = min(LCs.idx.neg2posRminusD(LCs.idx.neg2posRminusD > 0));
                LCs.idx.nextleastposRminusD = find(solved(i).TS.RminusD == LCs.idx.nextleastposRminusDval);
                % indexing the values for negative and positive R-D
                LCs.idx.nt = 0;
                LCs.idx.pt = 0;
                while LCs.idx.RminusDfit(solved(i).TS.t(LCs.idx.mostnegRminusD + LCs.idx.nt)) >= 0
                    LCs.idx.nt = LCs.idx.nt + 1;
                    if LCs.idx.mostnegRminusD + LCs.idx.nt > length(solved(i).TS.t)
                        solved(i).MC.tRminusD0 = NaN;
                        break
                    end
                end
                while LCs.idx.RminusDfit(solved(i).TS.t(LCs.idx.nextleastposRminusD + LCs.idx.pt)) <= 0
                    LCs.idx.pt = LCs.idx.pt + 1;
                    if LCs.idx.nextleastposRminusD + LCs.idx.pt > length(solved(i).TS.t)
                        solved(i).MC.tRminusD0 = NaN;
                        break
                    end
                end
            end
        end
        % find the zero time for R - D
        if ~isfield(solved(i).MC, 'tRminusD0')
            solved(i).MC.tRminusD0 = fzero(@(t) LCs.idx.RminusDfit(t), [solved(i).TS.t(LCs.idx.mostnegRminusD + LCs.idx.nt), solved(i).TS.t(LCs.idx.nextleastposRminusD + LCs.idx.pt)], optimset('Display','off'));
        end
    else
        solved(i).MC.tRminusD0 = NaN;
    end
    
    % conditional on direction of model
    if strcmp('rev', LCs.settings.timeDirection)
        % -- time [Gyr]
        solved(i).MC.t0 = solved(i).TS.t(end);
        
        % -- temperature [degC]
        solved(i).MC.T0 = solved(i).TS.T(end);
        solved(i).MC.Tpd = solved(i).TS.T(1);

        % -- upper mantle water concentration [ppm]
        solved(i).MC.Xum0 = solved(i).TS.Xum(end);
        solved(i).MC.Xumpd = solved(i).TS.Xum(1);
        
        % -- plate water concentration [ppm]
        solved(i).MC.Xppd = solved(i).TS.Xp(1);

        % -- regassing efficiency [%] 
        solved(i).MC.fRpd = solved(i).TS.fR(1);

        % -- plate thickness [km]
        solved(i).MC.deltapd = solved(i).TS.delta(1);

        % -- plate speed [cm/yr]
        solved(i).MC.upd = solved(i).TS.u(1);

        % -- total heat flow out [TW]
        solved(i).MC.Qspd = solved(i).TS.Qs(1);
        
        % -- calculated upper mantle viscosity [Pa s]
        solved(i).MC.eta_umpd = solved(i).TS.eta_um(1);

        % -- present-day surface ocean mass [OM]
        solved(i).MC.oceanmasspd = solved(i).TS.oceanmass(1);

        % -- present-day R - D
        solved(i).MC.RminusDpd = solved(i).TS.RminusD(1);

        % -- present-day Urey Ratio
        solved(i).MC.Ureypd = solved(i).TS.UreyRatio(1);
    else
        % -- temperature [degC]
        solved(i).MC.T0 = solved(i).TS.T(1);
        solved(i).MC.Tpd = solved(i).TS.T(end);

        % -- upper mantle water concentration [ppm]
        solved(i).MC.Xum0 = solved(i).TS.Xum(1);
        solved(i).MC.Xumpd = solved(i).TS.Xum(end);
	solved(i).MC.fH2Oumpd = solved(i).TS.fH2Oum(end);

        % -- plate water concentration [ppm]
        solved(i).MC.Xppd = solved(i).TS.Xp(end);

        % -- regassing efficiency [%] 
        solved(i).MC.fRpd = solved(i).TS.fR(end);
        
        % -- plate thickness [km]
        solved(i).MC.deltapd = solved(i).TS.delta(end);

        % -- plate speed [cm/yr]
        solved(i).MC.upd = solved(i).TS.u(end);

        % -- total heat flow out [TW]
        solved(i).MC.Qspd = solved(i).TS.Qs(end);
        
        % -- calculated upper mantle viscosity [Pa s]
        solved(i).MC.eta_umpd = solved(i).TS.eta_um(end);

        % -- present-day surface ocean mass [OM]
        solved(i).MC.oceanmasspd = solved(i).TS.oceanmass(end);

        % -- present-day R - D
        solved(i).MC.RminusDpd = solved(i).TS.RminusD(end);
        
        % -- present-day Urey Ratio
        solved(i).MC.Ureypd = solved(i).TS.UreyRatio(end);
    end
end
toc

%% MOVING DATA FROM solved TO ECs
% move all the MC variables over
ECs.solved.MCfieldnames = fieldnames(solved(1).MC);
ECs.solved.MC = [solved.MC];

%% CHECKING FOR HOW MANY MODELS FELL WITHIN UNCERTAINTIES
% initializing for successes
ECs.solved.successes = false(size(ECs.solved.MC));
% general stuff for T, OM, u, Qs, fR
if strcmp('rev', ECs.settings.timeDirection)
    myplots.solved.labelsBestOf = {'t', 'T', 'OM', 'u', 'Qs', 'fH2Oum'};
else
    myplots.solved.labelsBestOf = {'T', 'OM', 'u', 'Qs', 'fH2Oum'};
end
for i = 1:length(myplots.solved.labelsBestOf)
    myplots.solved.(['feasible', myplots.solved.labelsBestOf{i}]) = 0;
end

% finding those
ECs.solved.nFullSims = 0;
for j = 1:ECs.MC.iterations
    if solved(j).MC.fullSim || strcmp('rev', ECs.settings.timeDirection)
        ECs.solved.nFullSims = ECs.solved.nFullSims + 1;
        % t
        if strcmp('rev', ECs.settings.timeDirection) && solved(j).MC.t0 <= 0.75
            myplots.solved.feasibleT = myplots.solved.feasiblet + 1;
        end
        % T
        if (solved(j).MC.Tpd >= ECs.range.Tpd(1)-ECs.conversion.degC2K && solved(j).MC.Tpd <= ECs.range.Tpd(end)-ECs.conversion.degC2K)
            myplots.solved.feasibleT = myplots.solved.feasibleT + 1;
        end
        % OM
        if (solved(j).MC.oceanmasspd >= 0.9 && solved(j).MC.oceanmasspd <= 1.1)
            myplots.solved.feasibleOM = myplots.solved.feasibleOM + 1;
        end
        % u
        if (solved(j).MC.upd >= 4 && solved(j).MC.upd <= 5)
            myplots.solved.feasibleu = myplots.solved.feasibleu + 1;
        end
        % Qs
        if (solved(j).MC.Qspd >= 30 && solved(j).MC.Qspd <= 34)
            myplots.solved.feasibleQs = myplots.solved.feasibleQs + 1;
        end
	% f_H2O,um
	if (solved(j).MC.fH2Oumpd >= 0.35 && solved(j).MC.fH2Oumpd <= 63.3)
	    myplots.solved.feasiblefH2Oum = myplots.solved.feasiblefH2Oum + 1;
	end	

        % qualifying successes
	myplots.successConditions = (solved(j).MC.Tpd >= ECs.range.Tpd(1)-ECs.conversion.degC2K && solved(j).MC.Tpd <= ECs.range.Tpd(end)-ECs.conversion.degC2K) && ...
           			    (solved(j).MC.upd >= 4 && solved(j).MC.upd <= 5) && ...
           			    (solved(j).MC.Qspd >= 30 && solved(j).MC.Qspd <= 34) && ...
           			    (solved(j).MC.oceanmasspd >= 0.9 && solved(j).MC.oceanmasspd <= 1.1) && ...
           			    (solved(j).MC.fRpd <= max(ECs.range.fR)*100) && ...
				    (solved(j).MC.fH2Oumpd >= 0.35 && solved(j).MC.fH2Oumpd <= 63.3);
        if myplots.successConditions    
	    ECs.solved.successes(j) = true;
        end
    end
end

SectionBreak()
disp(['SEED AND GENERATOR USED: ', num2str(ECs.settings.seed.Seed), ', ', ECs.settings.seed.Type])
disp(['TOTAL COMPLETE REALIZATIONS: ', num2str(ECs.solved.nFullSims)])
for i = 1:length(myplots.solved.labelsBestOf)
    disp(['Fraction of complete realizations with feasible ',myplots.solved.labelsBestOf{i},': ', num2str(100 * myplots.solved.(['feasible',myplots.solved.labelsBestOf{i}]) / ECs.solved.nFullSims),'%'])
end
disp(['TOTAL SIMULTANEOUS FEASIBLE REALIZATIONS: ',num2str(sum(ECs.solved.successes))])

%% REMOVE REALIZATIONS THAT AREN'T FULL
if sum(ECs.solved.successes) > 0
    ECs.solved.filter = ~ECs.solved.successes; 
else
    ECs.solved.filter = [ECs.solved.MC(:).fullSim];
    ECs.solved.filter = ~ECs.solved.filter;
end
ECs.solved.successes(ECs.solved.filter) = [];
ECs.solved.MC(ECs.solved.filter) = [];

%% SAVE THE DATA
SectionBreak()

% should the data be saved at all?
if strcmp('base', ECs.settings.modelType) || sum(ECs.solved.successes) > 0
    ECs.settings.saveData = 'y';
else
    ECs.settings.saveData = 'n';
    disp('MODEL NOT SAVED...')
end

% what data to save and where to save it
if ECs.settings.saveData == 'y'
    % create an appropriate filepath...
    ECs.settings.filePath = fileparts(mfilename('fullpath'));
    
    % time direction
    if strcmp('fwd', ECs.settings.timeDirection)
        ECs.settings.filePath = [ECs.settings.filePath,'/models/fwd/'];
    else
        ECs.settings.filePath = [ECs.settings.filePath,'/models/rev/'];
    end
    
    % model type
    if strcmp('base', ECs.settings.modelType)
        ECs.settings.filePath = [ECs.settings.filePath,'base/'];
    else
        ECs.settings.filePath = [ECs.settings.filePath,'var/'];
    end
    
    % create the folder if it doesn't exist
    if ~exist(ECs.settings.filePath,'dir')
        mkdir(ECs.settings.filePath)
    end

    % if this is not a baseline model...
    ECs.settings.varFileName = '';
    if strcmp('var', ECs.settings.modelType)
        % variable fR
        if ECs.settings.varyfR == 'y'
            ECs.settings.varFileName = [ECs.settings.varFileName, 'VarfR_'];
        end

        % variable phiRum
        if ECs.settings.varyphiRum == 'y'
            ECs.settings.varFileName = [ECs.settings.varFileName, 'VarphiRum_'];
        end
    end

    % create the file name
    ECs.settings.fileName = [upper(ECs.settings.timeDirection), '_', ...
                             upper(ECs.settings.modelType), '_', ...                         
                             upper(ECs.settings.viscDep), '_', ...
                             upper(ECs.settings.viscWaterStrength), '_', ...
                             ECs.settings.varFileName, ...
                             'MC'];
    
    % saving MC data
    disp('SAVING MC DATA...')
    save([ECs.settings.filePath, ECs.settings.fileName], 'ECs')

    % only save the timeseries data for the baseline cases
    if strcmp('base', ECs.settings.modelType) || ECs.settings.varyfR == 'y'
        ECs.solved = rmfield(ECs.solved, 'MC');
        ECs.solved.TS = [solved.TS];
        ECs.solved.TS(ECs.solved.filter) = [];
        disp('SAVING TS DATA...')
        save([ECs.settings.filePath, [upper(ECs.settings.timeDirection), '_', ...
                                      upper(ECs.settings.modelType), '_', ...
                                      upper(ECs.settings.viscDep), '_', ...
                                      upper(ECs.settings.viscWaterStrength), '_', ...
                                      ECs.settings.varFileName, ...
                                      'TS']], 'ECs', '-v7.3')
        
    end
    
    disp('SUCCESSFULLY SAVED!')
    SectionBreak()
end


%% SECTIONBREAK()
function SectionBreak()
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
end

%% CONSERVATION EQUATIONS
function dYdt = ConservationEquations(t, Y, lcs)
    % Y(1) = T(t)
    % Y(2) = X_{um}(t)
    
    % MELT DEPTH
    z = lcs.D.z1*Y(1) + lcs.D.z2*Y(2) + lcs.D.z3;

    % HEAT PRODUCTION
    H = lcs.avg.rho * lcs.avg.V_M * ...
       (lcs.params.C_U238 * lcs.consE.H_U238 * exp(log(2) * (lcs.pd.ageEarth - t) / lcs.consE.tau_U238) + ...
        lcs.params.C_U235 * lcs.consE.H_U235 * exp(log(2) * (lcs.pd.ageEarth - t) / lcs.consE.tau_U235) + ...
        lcs.params.C_Th232 * lcs.consE.H_Th232 * exp(log(2) * (lcs.pd.ageEarth - t) / lcs.consE.tau_Th232) + ...
        lcs.params.C_K40 * lcs.consE.H_K40 * exp(log(2) * (lcs.pd.ageEarth - t) / lcs.consE.tau_K40));
    
    % VISCOSITY
    if strcmp('conc', lcs.settings.viscDep)
        H2Oconc = Y(2) / lcs.visc.Xc;
        eta_um = lcs.scale.eta0 * (((H2Oconc)^(-lcs.visc.r)) * exp(lcs.params.E / (lcs.avg.Rg * Y(1))));
    else
        H2Ofug = exp(lcs.visc.fugC0 + ...
                     lcs.visc.fugC1 * log(lcs.visc.ppm2COH * Y(2) / lcs.conversion.ppm2massfrac) + ...
                     lcs.visc.fugC2 * (log(lcs.visc.ppm2COH * Y(2) / lcs.conversion.ppm2massfrac))^2 + ...
                     lcs.visc.fugC3 * (log(lcs.visc.ppm2COH * Y(2) / lcs.conversion.ppm2massfrac))^3);
        eta_um = lcs.scale.eta0 * (((H2Ofug)^(-lcs.visc.r)) * exp(lcs.params.E / (lcs.avg.Rg * Y(1))));
    end
    
    % RAYLEIGH NUMBER
    Ra = ((lcs.avg.alpha * lcs.avg.rho * lcs.avg.g) * (Y(1) - lcs.pd.Tsurf) * lcs.avg.R_Em^3) / (lcs.avg.kappa * eta_um);

    % PLATE SPEED
    u = lcs.scale.u * (lcs.avg.kappa / lcs.avg.R_Em) * (lcs.avg.L / (pi * lcs.avg.R_Em))^(lcs.params.beta) * Ra^(2*lcs.params.beta);

    % PLATE THICKNESS
    delta = 2 * sqrt((lcs.avg.kappa * lcs.avg.L) / u);

    % TOTAL HEAT FLOW THROUGH SURFACE AREA
    Qs = 2 * lcs.scale.Qs * lcs.avg.A_O * lcs.avg.k * (Y(1) - lcs.pd.Tsurf) * sqrt(u / (pi * lcs.avg.L * lcs.avg.kappa));

    % DEGASSING EXCHANGE RATE
    D = lcs.avg.A_O * lcs.D.fD * (z / lcs.avg.L) * u * lcs.avg.rhoum * Y(2);

    % PLATE WATER CONCENTRATION
    Xp = lcs.params.Xp;

    % REGASSING EXCHANGE RATE
    % -- when we're considering a variable fR
    if lcs.settings.varyfR == 'y'
        % need u in [cm/yr], a in [Myr], and T in [degC]
        a = ((delta / lcs.params.W(5))^2) / lcs.avg.kappa * lcs.conversion.sec2yr * 1e-6;
        W = 1e5 * (lcs.params.W(1) * (u * lcs.conversion.yr2sec * 100) + ...
                   lcs.params.W(2) * a + ...
                   lcs.params.W(3) * (Y(1) - lcs.conversion.degC2K) + ...
                   lcs.params.W(4)*lcs.scale.MagniCScale);
        W0 = Xp * lcs.avg.rhop * (lcs.avg.A_O / lcs.avg.T);
        if (lcs.scale.fR0 * (W / W0)) >= 1
            fR = 1-eps;
        elseif (lcs.scale.fR0 * (W / W0)) <= 0
            fR = eps;
        else
            fR = (lcs.scale.fR0 * (W / W0));
        end
    else
        fR = lcs.params.fR;
    end
    R = lcs.avg.A_O * fR * (delta / lcs.avg.L) * u * lcs.avg.rhop * Xp;

    % ENERGY CONSERVATION
    dTdt = (H - Qs) / (lcs.avg.c_p * lcs.avg.rho * lcs.avg.V_M);

    % MASS CONSERVATION
    if lcs.settings.varyphiRum == 'y'
        dXumdt = (R*lcs.params.phiRum - D) / (lcs.avg.rhoum * lcs.avg.Vum);
    else
        dXumdt = (R - D) / (lcs.avg.rhoum * lcs.avg.Vum) * lcs.params.fum;
    end

    % dY/dt
    dYdt = [dTdt; dXumdt];
end

