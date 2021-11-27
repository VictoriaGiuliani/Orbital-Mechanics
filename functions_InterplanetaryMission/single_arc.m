function [dV1, dV2, V_m, V_p] = single_arc (t1,t2,id1,id2)
% SINGLE HELIOCENTRIC TRANSFER LEG: calculates the heliocentric orbital 
%   arc going from one planet to another planet in a set time.
%   In other words it's a Lambert arc designer where the initial and final
%   positions are taken from the ephemerides of the planets at departure
%   and arrival times
%__________________________________________________________________________   
% PROTOTYPE:
%    [dV1, dV2, V_m, V_p] = single_arc (t1,t2,id1,id2)
% 
% INPUT:
%   t_1[1]    time in MJD2000, for the departure of the s/c from the
%             first specified planet.                               [days]
%   t_2[1]    time in MJD2000, for the arrival of the s/c to the
%             second specified planet.                              [days]
%   id1[1]    departure planet identifier                           [-]
%   id2[1]    arrival planet identifier                             [-]
%
% OUTPUT:
%   dV1[1]    DeltaV given by the difference between the first 
%             planet's heliocentric velocity and the starting velocity
%             of the Lambert arc                                    [km/s]
%   dV2[1]    DeltaV given by the difference between the second 
%             planet's heliocentric velocity and the end velocity
%             of the Lambert arc                                    [km/s]
%   V_m[1X3]  Heliocentric end velocity of the Lambert arc          
%             (at the second planet)                                [km/s]
%   V_p[1X3]  Heliocentric starting velocity of the Lambert arc     
%             (at the first planet)                                 [km/s]
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find first planet's position and velocity from ephemerides
[kep1,muSun] = uplanet(t1, id1);
[r1,v1] = kep2car(kep1, muSun);

% find second planet's position and velocity from ephemerides
[kep2,muSun] = uplanet(t2, id2);
[r2,v2] = kep2car(kep2, muSun);

% calculate the transfer time from first to second planet in seconds
Dt = (t2 - t1)*24*3600;

% calculate only physically feasible transfers (no negative transfer times)
if Dt > 0
    % find initial and final velocities of the Lambert arc
    [~, ~, ~, ~, V_p, V_m, ~, ~] = lambertMR( r1, r2, Dt, muSun, 0, 0, 0, 2 );
    % find the two DeltaVs
    dV1 = norm(v1 - V_p);
    dV2 = norm(v2 - V_m);
else
    % if the transfer isn't physically feasible, output NaN
    V_m = NaN;
    V_p = NaN;
    dV1 = NaN;
    dV2 = NaN;
end
