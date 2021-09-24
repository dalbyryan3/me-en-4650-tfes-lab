function [rho,Cp] = WaterProperties(T)
%--------------------------------------------------------------------
% function [rho,Cp] = WaterProperties(T)
%
% Properties of saturated liquid water at atmospheric pressure. Useful
% in analyzing heat exchanger lab data.
%
% INPUTS:
%    T      Temperature in oC
%
% OUTPUTS:
%    rho    density in kg/m^3
%    Cp     specific heat in J/(kg*K)
%
% M Metzger
%--------------------------------------------------------------------

% initialize output variables
rho=zeros(size(T));
Cp=zeros(size(T));

% Temperature Data for Specific Volume (oC)
T_vec_V = [0.01,10:10:360];   

% Specific Volume (m^3/kg)
V_vec = [1,1,1.002,1.004,1.008,1.012,1.017,1.023,1.029,1.036,1.043,...
    1.052,1.06,1.07,1.08,1.091,1.102,1.114,1.127,1.141,1.157,1.173,...
    1.19,1.209,1.229,1.252,1.276,1.303,1.333,1.366,1.404,1.447,...
    1.499,1.56,1.638,1.741,1.895]*1e-3;

% Temperature Data for Specific Heat (oC)
T_vec_Cp = [0.01,10,20,25,30:10:120,140:20:360];

% Specific heat (J/kg/K)
Cp_vec = [4.2199,4.1955,4.1844,4.1816,4.1801,4.1796,4.1815,4.1851,...
    4.1902,4.1969,4.2053,4.2157,4.2283,4.2435,4.2826,4.3354,4.4050,...
    4.4958,4.6146,4.7719,4.9856,5.2889,5.7504,6.5373,8.2080,15.004]*1e3;

% interpolate to find properties at each given temperature
for (j=1:length(T))
    V = interp1(T_vec_V,V_vec,T(j));        %specific volume
    Cp(j) = interp1(T_vec_Cp,Cp_vec,T(j));  %specific heat
    rho(j) = 1/V;                           %density
end

