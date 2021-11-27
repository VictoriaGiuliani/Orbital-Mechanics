function [i, OM, om] = findplane(u,vinfM,delta)
% FINDPLANE function that computes the 3D orientation of the hyperbolic
% orbit used in a flyby
% _________________________________________________________________________  
% PROTOTYPE:
%   [i, OM, om] = findplane(u,vinfM,delta)
%
% INPUT:
%   u       [3x1]   Direction around which the velocity on the orbit
%                   rotates, it's perpendicular to the orbital plane and
%                   points in the direction indicated by the right hand
%                   rule.
%                   It's a versor.                                      [-]
%   vinfM   [3x1]   Incoming velocity of the s/c at time that tends to
%                   minus infinity. Given in the planet-centered reference
%                   frame.                                           [km/s]
%   delta   [1]     Deflection/Turn angle of the hyperbola.           [rad]
%
% OUTPUT:
%   i       [1]     Inclination of the orbital plane                  [rad]
%   OM      [1]     RAAN of the orbital plane                         [rad]
%   om      [1]     Anomaly of the pericenter of the orbit            [rad]
% _________________________________________________________________________
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi
% _________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Finding inclination
i = acos(u(3));

% Finding the versor N that points to the ascending node
k = [0 0 1]';
N = cross(k, u);
N = N / norm(N);
if isnan(N)
    N = [1; 0; 0];
end

% Finding RAAN as angle between N and the vernal equinox line
if N(2) >= 0
	OM = acos(N(1));
elseif N(2) < 0
	OM = 2 * pi - acos(N(1));
end

% Finding the angle beta of the hyperbola
beta = (pi-delta)/2;

% Finding the direction of the eccentricity vector of the hyperbola. It's
% the direction of vinfM rotated by beta using Rodrigues' rotation formula.
ee = -(vinfM * cos(beta) + cross(u,vinfM) * sin(beta)...
        + u * dot(u,vinfM) * (1-cos(beta)));
e = norm(ee); % Note: this is NOT the eccentricity of the hyperbola.

% Finding anomaly of the perigee as angle between eccentricity and nodes
% directions.
if ee(3) >= 0
	om = acos(dot(N, ee) / e);
elseif ee(3) < 0
	om = 2 * pi - acos(dot(N, ee) / e);
end

end
