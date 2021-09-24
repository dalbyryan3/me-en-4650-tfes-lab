function [rho,mu,k,Cp] = AirProperties(T,P)
%--------------------------------------------------------------------
% function [rho,nu,k,Cp] = AirProperties(T,P)
%
% Properties of dry air based on the ideal gas law. Note, for ideal
% gases, Cp, k, and mu are independent of pressure.
%
% INPUTS:
%    T      Temperature in Kelvin
%    P      Atmospheric pressure in Pascal
%
% OUTPUTS:
%    rho    density in kg/m^3
%    mu     absolute viscosity in kg/(m*s)
%    k      thermal conductivity in W/(m*K)
%    Cp     specific heat in J/(kg*K)
%
% M Metzger
%--------------------------------------------------------------------

% gas constant for dry air
Rair = 287.05;  %J/kg*K

% density based on the ideal gas law (kg/m^3)
rho = P/(Rair*T);

% Temperature Data (K)
T_vec = [-150,-100,-50:10:-10,0:5:50,60:10:100,120:20:200,...
    250:50:500,600:100:1000,1500,2000]+273.15;   

% absolute viscosity (kg/m*s)
mu_vec = [0.8636,1.189,1.474,1.527,1.579,1.630,1.680,...
    1.729,1.754,1.778,1.802,1.825,1.849,1.872,1.895,1.918,1.941,...
    1.963,2.008,2.052,2.096,2.139,2.181,2.264,2.345,2.420,2.504,2.577,...
    2.76,2.934,3.101,3.261,3.415,3.563,3.846,4.111,4.362,4.6,...
    4.826,5.817,6.630]*1e-5;

% Thermal conductivity (W/m/K)
k_vec = [0.01171,0.01582,0.01979,0.02057,0.02134,0.02211,0.02288,...
    0.02364,0.02401,0.02439,0.02476,0.02514,0.02551,0.02588,0.02625,...
    0.02662,0.02699,0.02735,0.02808,0.02881,0.02953,0.03024,0.03095,...
    0.03235,0.03374,0.03511,0.03646,0.03779,0.04104,0.04418,0.04721,...
    0.05015,0.05298,0.05572,0.06093,0.06581,0.07037,0.07465,...
    0.07868,0.09599,0.11113];

% Specific heat (J/kg/K)
Cp_vec = [983,966,999,1002,1004,1005,1006,...
    1006,1006,1006,1007,1007,1007,1007,1007,1007,1007,1007,1007,...
    1007,1008,1008,1009,1011,1013,1016,1019,1023,1033,1044,1056,...
    1069,1081,1093,1115,1135,1153,1169,1184,1234,1264];

% interpolate to find properties at given temperature
mu = interp1(T_vec,mu_vec,T);
k = interp1(T_vec,k_vec,T);
Cp = interp1(T_vec,Cp_vec,T);
