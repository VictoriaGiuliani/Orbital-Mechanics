function [rp, dv_GA, delta_m, delta_p, a_m, a_p, e_m, e_p] = GApow(V_m, V_p, idGA, tGA)
% POWERED GRAVITY ASSIST: from the incoming and outgoing heliocentric
%   velocities of the flyby planet, it calculates the orbits inside the
%   SOI of the planet.
%   The impulsive manoeuvre is assumed at the hyperbola's pericentre and 
%   in the direction of flight. Therefore the flyby has two hyperbolas,
%   patched at their pericentre.
%   The function checks if the s/c will impact the planet.
%   NOTE: heliocentric velocites will use capital V, while planet-centered
%   velocities will use small v
%__________________________________________________________________________   
% PROTOTYPE:
%    [rp, dv_GA] = GApow(V_m, V_p, idGA, tGA)
% or [rp,dv_GA,delta_m,delta_p,a_m,a_p,e_m,e_p] = GApow(V_m, V_p, idGA, tGA)
%
% INPUT:
%   V_m[1X3]  heliocentric velocity of the incoming s/c at infinity [km/s]
%   V_p[1X3]  heliocentric velocity of the outgoing s/c at infinity [km/s]
%   idGA[1]   flyby planet identifier                               [-]
%   tGA[1]    MJD2000 time at which the flyby happens
%             (flyby assumed as instantaneous from the heliocentric
%             point of view)                                        [days]
%
% OUTPUT:
%   rp[1]     Radius of the pericentre of the two hyperbolas, if
%             the flyby is unfeasible, it's NaN                     [km]
%   dv_GA[1]  DeltaV required for the powered part of the gravity
%             assist                                                [km/s]
%   delta_m[1]Deflection/turn angle of the incoming hyperbolic arc  [rad]
%   delta_p[1]Deflection/turn angle of the outgoing hyperbolic arc  [rad]
%   a_m[1]    Semi-major axis of the incoming hyperbolic arc        [km]
%   a_p[1]    Semi-major axis of the outgoing hyperbolic arc        [km]
%   e_m[1]    Eccentricity of the incoming hyperbolic arc           [-]
%   e_p[1]    Eccentricity of the outgoing hyperbolic arc           [-]
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% atmosphere thickness of the planet, needs to be specified by the user in
% km
global h_atm 

% calculate parameters of the gravity assist planet
muGA = astroConstants(idGA + 10);
R_GA = astroConstants(20 + idGA);
% minimum radius of perientre to achieve physical feasiblilty
rp_min = R_GA + h_atm;
if isempty(rp_min)
    rp_min = R_GA;
end

% Find velocity of the flyby planet (assumed constant during flyby)
[kep, muS] = uplanet(tGA, idGA);
[~,V] = kep2car(kep, muS); 

% Find arrival and departure velocities at infinite wrt the planet
v_inf_m = V_m - V;
v_inf_p = V_p - V;

% Find the magnitude of these velocities
v_inf_m_n = norm(v_inf_m);
v_inf_p_n = norm(v_inf_p);

% Find the turning angle of the total powered transfer
delta = acos(dot(v_inf_m, v_inf_p)/(v_inf_m_n*v_inf_p_n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to solve the upcoming implicit function used to calculate rp, 
% an estimate of rp is necessary:

% find the semi major axes of the two hyperbolas
a_m = - muGA / (v_inf_m_n)^2;
a_p = - muGA / (v_inf_p_n)^2;

% find two estimates for the Impact Parameters of the inital and final
% hyperbolas, the total delta is assumed for both hyperbolas, but this is
% just an approximation
Delta1_guess = - a_m / tan(delta/2);
Delta2_guess = - a_p / tan(delta/2);

% Using the approximations of impact parameters, guesses for the pericentre
% radii of each hyperbola are found
rp_guess1 = (-muGA + sqrt(muGA^2 + Delta1_guess^2*v_inf_m_n^4))/v_inf_m_n^2;
rp_guess2 = (-muGA + sqrt(muGA^2 + Delta2_guess^2*v_inf_p_n^4))/v_inf_p_n^2;

% The initial guess of the implicit function will be the average of the two
% radii
rp_guess  = 0.5*(rp_guess1+rp_guess2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check before solving the implicit function whether the s/c will impact
% the planet. 
% For a rougher but quicker and more stable solution substitute || with &&.
% With && some borderline cases may not be computed (this doesn't influence
% the local minima search).
if rp_guess1 > rp_min && rp_guess2 > rp_min
    
    %Computation of rp through an implicit function
    e_m = @(rp) 1 + rp*(v_inf_m_n)^2/muGA;
    e_p = @(rp) 1 + rp*(v_inf_p_n)^2/muGA;
    f = @(rp) asin(1./e_m(rp)) + asin(1./e_p(rp)) - delta; 
    options = optimset('TolFun',1e-12);
    rp = fzero(f, rp_guess,options);
    
    
    % Compute the velocities of the two hyperbolic arcs at pericentre 
    % and the required powered DeltaV at pericentre
    vp_m = sqrt(2*(0.5*(norm(v_inf_m))^2 + muGA/rp));
    vp_p = sqrt(2*(0.5*(norm(v_inf_p))^2 + muGA/rp));

    dv_GA = abs(vp_m - vp_p);   
    
    % Recheck feasibility with real value of rp
    if rp < rp_min
        % What happens if transfer is unfeasible
        dv_GA = NaN;
        rp = NaN;
    end

    % Characteristics of the two hyperbolic arcs, semimajor axes will be
    % equal to the ones calculated at lines 74-75
    a_m = rp/(1 - e_m(rp));
    a_p = rp/(1 - e_p(rp));

    e_m = e_m(rp);
    e_p = e_p(rp);
    
    delta_m = 2*asin(1/e_m);
    delta_p = 2*asin(1/e_p);
else
    % What happens if transfer is unfeasible
    dv_GA = NaN;
    rp = NaN;
    a_m=NaN;a_p=NaN;e_m=NaN;e_p=NaN;delta_m=NaN;delta_p=NaN;
end
end


















