function [r, v] = kep2car(kep, mu)
% KEPLERIAN TO CARTESIAN: % Conversion from Keplerian elements to 
%   Cartesian coordinates
%__________________________________________________________________________   
% PROTOTYPE:
%    [r, v] = kep2car(kep, mu)
% 
% INPUT:
%   kep   [1x6]    Vector containing the following parameters:
%         a     [1x1]    Semi-major axis            [km]
%         e     [1x1]    Eccentricity               [-]
%         i     [1x1]    Inclination                [rad]
%         OM    [1x1]    RAAN                       [rad]
%         om    [1x1]    Pericentre anomaly         [rad]
%         th    [1x1]    True anomaly               [rad]
%         mu    [1x1]    Gravitational parameter    [km^3/s^2]
% 
% OUTPUT:
%   r     [1x3]    Position vector            [km]
%   v     [1x3]    Velocity vector            [km/s]
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a  = kep(1);
e  = kep(2);
i  = kep(3);
OM = kep(4);
om = kep(5);
th = kep(6);

p = a * (1 - e^2);                     % Semi-latus rectum
int_r = p / (1 + e*cos(th));           % r intensity

r_pf = int_r * [cos(th) sin(th) 0];              % r in perifocal coordinate system
v_pf = sqrt(mu/p) * [-sin(th) e+cos(th) 0];      % v in perifocal coordinate system

R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];     % Matrix for rotation around K axis of OM rad 
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];          % Matrix for rotation around I axis of i rad 
R3_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];     % Matrix for rotation around K axis of om rad 

% Change of coordinates matrix: from PF to ECI
R = R3_OM' * R1_i' * R3_om';
vv_1 = sqrt(mu / p) * [- sin(th); e + cos(th); zeros(size(th))];

r = (R * r_pf')';                  % r in Cartesian coordinates
v = (R * v_pf')';                  % v in Cartesian coordinates
