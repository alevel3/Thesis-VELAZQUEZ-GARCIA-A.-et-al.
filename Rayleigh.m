function [Ext_air]=Rayleigh(Pression,Temperature,lambda,CCO2)
% BODHAINE et al. 1999 Journal of atmospheric and oceanic technology 
% Exemple : 
% CCO2 = 360; % ppm
% Pression = 1013.25;
% Temperature= 288.15;
lambda = lambda./1000; % um
CCO2 = CCO2/10000; % units should be first in ppm ex: CCO2 = 360ppm

if Temperature <200
    Temperature = Temperature +273.15;
end


%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Avogadro's number
Na=6.02214*1e23;        %%unit: 1/mol
%% Gas constant
Ra=8.314472;            %%unit: J/K/mol
%% Boltzmann's constant
kB=1.3806504*1e-23;     %%unit: J/K
%% Standard pressure and temperature
Ps=1013.25;             %%unit: hPa
Ts=288.15;              %%unit: K
%% 
m_air = 15.0556*CCO2/100 + 28.9595; %% gm mol-1
g = 9.81 * 100; % (cm s-2)

%% refractive index : 
% We recommend starting with Peck and Reeder?s
% (1972) formula for the refractive index of dry air with
% 300 ppm CO2 concentration:
n_air_300 =  (8060.51 + 2480990./(132.274 - lambda^-2) + 17455.7./(39.32957 -lambda^-2))/1e8 +1; 
n_air_CCO2 = ((1+0.54*(CCO2-0.0003))*(n_air_300-1))+1;


Ns=(Na*Pression*1e2)/(Ra*Temperature); %%//unit: m-3
FN2 = 1.034 +3.17*1e-4./lambda^2;
FO2 = 1.096 + 1.385*10^-3./lambda^2 + 1.448*10^-4./lambda^4;
FAr = 1.00;
FCO2 = 1.15;
Fair = (78.084*FN2 + 20.946*FO2 +0.934*FAr + CCO2*FCO2)./(78.084 + 20.946 +0.934 + CCO2);

Sigma_air = 24 * pi^3*(n_air_CCO2^2-1)^2*(6+3*Fair)./(lambda^4*Ns^2*(n_air_CCO2^2+2)^2*(6-7*Fair));

Sigma_air_test = 1e-28*(1.0455996 - 341.2961*lambda^-2 - 0.90230850*lambda^2)./(1+ 0.0027059889*lambda^-2 - 85.968563*lambda^2);

Ext_air = Sigma_air_test * Pression*1000 * Na./(m_air *g) .*100; %Mm-1
end