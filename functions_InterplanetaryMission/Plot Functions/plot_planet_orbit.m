function plot_planet_orbit(t_start,t_end,id)
% PLOT PLANET ORBIT: plots the 3D orbit of the considered planet on the
%   currently open figure, in the ecliptic frame of reference.
%   Highlits the part of the orbit spent between t_start and t_end.
%__________________________________________________________________________   
% PROTOTYPE:
%    plot_planet_orbit(t_start,t_end,id)
% 
% INPUT:
%   t_start[1]  time in MJD2000, specifies where the highlited part
%               of the orbit should start                           [days]
%   t_end[1]    time in MJD2000, specifies where the highlited part
%               of the orbit should end. If it's equal to t_start,
%               then no highlited part will be plotted              [days]
%   id[1]       planet identifier                                   [-]
%
% OUTPUT:
%   figure containing the 3D plot
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The plot is in astronomical units:
AU = astroConstants(2);

% Find inital position and velocity of the planet with ephemerides
[kep,mu_S]= uplanet(t_start, id);
[r0,v0] = kep2car( kep, mu_S );

% Instead of using ephemeris at each point, the planet's orbit can just be
% integrated, since just one period is considered

% Integrate and plot the highlited part of the orbit
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
if t_end>t_start
    tspan = linspace(t_start,t_end,1000)*24*3600;
    [~,y] = ode113(@(t,y) ODE_two_body(t,y,mu_S),tspan, [r0';v0'], options);
    plot3(y(:,1)/astroConstants(2),y(:,2)/astroConstants(2),y(:,3)/astroConstants(2),'LineWidth',2)
    hold on
end

% Integrate and plot the whole orbit
t_period = 2*pi*sqrt(kep(1)^3 / mu_S) /24/3600;
tspan = linspace(t_start,t_start + t_period,1000)*24*3600;
[~,y] = ode113(@(t,y) ODE_two_body(t,y,mu_S),tspan, [r0';v0'], options);
plot3(y(:,1)/AU,y(:,2)/AU,y(:,3)/AU,'k--','LineWidth',1)


end