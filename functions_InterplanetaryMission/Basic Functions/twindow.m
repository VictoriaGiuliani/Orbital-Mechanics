function [t1,t2] = twindow(id1,id2,weight_min,weight_max)
% TIME WINDOW: Find the hohman transfer time in the 2D circular case 
%   between planets identified by the inputs.
%   Then estimate the time window used for the optimum search in the
%   interplanetary transfer. The time window's extremes are scaled form the
%   Hohmann transfer time with either user set weights or standard 0.1 and 
%   2 weights taken from: https://bit.ly/386VLgk   
%__________________________________________________________________________   
% PROTOTYPES:
%    [t1,t2] = twindow(id1,id2,weight_min,weight_max)
% or [t1,t2] = twindow(id1,id2)
% 
% INPUT:
%   id1[1]          departure planet identifier                     [-]
%   id2[1]          arrival planet identifier                       [-]
%   weight_min[1]   minimum time limit with respect to the Hohmann
%                   transfer time                                   [-]
%   weight_max[1]   maximum time limit with respect to the Hohmann
%                   transfer time                                   [-]
%
% OUTPUT:
%   t1[1]           Lower search region boundary in days            [days]  
%   t2[1]           Upper search region boundary in days            [days]
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Take keplerian elements of planets from ephemerides, since this is a
% rough estimate time 0 is chosen
kep1 = uplanet(0,id1);
kep2 = uplanet(0,id2);
% Gravitational parameter of the Sun
muS  = astroConstants(4);

% radii of the assumed circular orbits
r1 = kep1(1); % supposed equal to a because of circular orbit hypothesis
r2 = kep2(1); % supposed equal to a a because of circular orbit hypothesis

% semi-major axis of the Hohmann transfer ellipse
aT = (r1 + r2)/2;

% Hohmann transfer time in days
t_hoh = pi * sqrt(aT^3 / muS);
t_hoh = t_hoh / (3600*24);

% Find outputs
if nargin == 2
    t1 = t_hoh * 0.1;
    t2 = t_hoh * 2.0;
elseif nargin == 4
    t1 = t_hoh * weight_min;
    t2 = t_hoh * weight_max;
else
    error([],'Wrong number of inputs, insert either 2 or 4 inputs.')
end


