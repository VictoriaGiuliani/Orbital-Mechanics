function kep = car2kep (r, v, mu)

%  CARTESIAN TO KEPLERIAN: Conversion from Cartesian coordinates to 
%    Keplerian elements 
%__________________________________________________________________________   
% PROTOTYPE:
%    kep = car2kep (r, v, mu)
% 
% INPUT:
%   r     [1x3]    Position vector            [km]
%   v     [1x3]    Velocity vector            [km/s]
%
% OUTPUT:
%   kep   [1x6]    Vector containing the following parameters:
%         a     [1x1]    Semi-major axis            [km]
%         e     [1x1]    Eccentricity               [-]
%         i     [1x1]    Inclination                [rad]
%         OM    [1x1]    RAAN                       [rad]
%         om    [1x1]    Pericentre anomaly         [rad]
%         th    [1x1]    True anomaly               [rad]
%         mu    [1x1]    Gravitational parameter    [km^3/s^2]
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    mu = 398600.433;                      
end

int_r = norm(r);     % r intensity
int_v = norm(v);     % v intensity

h = cross(r,v);      % Specific angolar momentum
int_h = norm(h);     % h intensity

i = acos( h(3) / int_h);   % Inclination = acos(hz/h)

e = 1/mu * ( ((int_v)^2 - mu/int_r) * r - dot(r,v) * v );   % Eccentricity vector
int_e = norm(e);                                            % e intensity

eps = 0.5 * (int_v)^2 - mu/int_r;       % Specific energy
a = - mu / (2*eps);                     % Semi-major axis

K = [0 0 1];             % K axis
N = cross(K,h);          % Nodal axis
int_N = norm(N);         % N intensity

% Definition of acos
if N(2)<0
        OM = 2*pi - acos( N(1) / int_N );    % RAAN if Ny < 0
else
        OM = acos( N(1) / int_N );           % RAAN if Ny >= 0
end

% Definition of acos
if e(3)<0
        om = 2*pi - acos( dot(N,e) / (int_N * int_e)  );    % Pericentre anomaly if ez < 0
else
        om = acos( dot(N,e) / (int_N * int_e));             % Pericentre anomaly if ez >= 0
end

% Radial velocity
% if vr > 0 => we are getting away from pericenter
% if vr < 0 => we are approaching pericenter
vr = dot(r,v) / int_r;    

% Definition of acos
if vr < 0
        th = 2*pi - acos( dot(e,r) / (int_r * int_e)  );    % True anomaly if vr < 0
else
        th = acos( dot(e,r) / (int_r * int_e));             % True anomaly if vr > 0
end

kep = [a int_e i OM om th];
