function [dV_tot,dV_dep,dV_arr,dV_ga] = GA_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3)
% INTERPLANETARY TRANSFER WITH GRAVITY ASSIST: calculates the
%   interplanetary transfer between a departure, gravity assist (which can 
%   be powered), and arrival planet, using a patched conic method.
%   No departure nor arrival hyperbolic trajectories are considered,
%   all manoeuvers (including the gravity assist) are assumed as
%   instantaneous.
%__________________________________________________________________________   
% PROTOTYPE:
%    [dV_tot,dV_dep,dV_arr,dV_ga] = GA_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3)
% 
% INPUT:
%   t_dep[1]  time in MJD2000, for the departure of the s/c.        [days]
%   t_ga[1]   time in MJD2000, for the gravity assist of the s/c,
%             it's assumed as instantaneous                         [days]
%   t_arr[1]  time in MJD2000, for the arrival of the s/c.          [days]
%   id1[1]    departure planet identifier                           [-]
%   id2[1]    gravity assist planet identifier                      [-]
%   id3[1]    arrival planet identifier                             [-]
%
% OUTPUT:
%   dV_tot[1] total DeltaV required for the complete transfer       [km/s]
%   dV_dep[1] DeltaV required for the departure from the first
%             planet, the departure orbits inside the SOI      
%             are neglected                                         [km/s]
%   dV_arr[1] DeltaV required for the arrival to the last
%             planet, the arrival orbits inside the SOI 
%             are neglected                                         [km/s]
%   dV_ga [1] DeltaV used in the powered gravity assist/flyby
%             of the second planet                                  [km/s]
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the departure DeltaV and the helicentric arrival velocity at
% the second planet (V_m)
[dV_dep, ~, V_m, ~] = single_arc(t_dep,t_ga,id1,id2); 

% calculate the arrival DeltaV and the helicentric departure velocity from
% the second planet (V_p)
[~, dV_arr, ~, V_p] = single_arc(t_ga,t_arr,id2,id3);

% calculate the DeltaV required from the s/c during the powered flyby, note
% that this doesn't include the DeltaV given by the gravity assist itself
% if the flyby is physically unfeasible, then dV_ga will be output as NaN
[~, dV_ga] = GApow(V_m, V_p, id2, t_ga); 

% sum up the DeltaVs
dV_tot = dV_dep + dV_ga + dV_arr;

end
