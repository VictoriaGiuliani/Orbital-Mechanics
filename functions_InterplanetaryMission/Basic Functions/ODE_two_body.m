function dy = ODE_two_body(~, y, mu)
% ODE TWO BODY PROBLEM: gives the derivatives (velocity and acceleration)
%   related to the state space formulation of the two body problem. 
%__________________________________________________________________________   
% PROTOTYPE:
%    dy = ODE_two_body(time,y,mu)
% 
% INPUT:
%   y[6x1]      position and velocity in cartesian coordinates     [km km km km/s km/s km/s]
%   mu[1]       planetary gravity constant                        [km^3/s^2]
%   
% OUTPUT:
%   velocity and acceleration as function of position and velocity
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=[ y(4:6); ...
    -mu./(norm(y(1:3)))^3.*y(1:3)];

end