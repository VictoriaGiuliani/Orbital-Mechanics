function plot_flyby(t_dep, t_ga, t_arr, id1, id2, id3)
% PLOT OF FLYBY: plots the powered flyby of the assigned 3 planet transfer 
%   (with the flyby being at the second planet).
%   The 3D plot will be centered at the planet, but the versors will be 
%   in the ecliptic reference frame.
%   The orbits are found under the same assumptions of the GApow function.
%   Optionally, the velocity triangles can be plotted.
%   Optionally, useful parameters can be printed on the command window.
%   To do so, uncomment the end of the function.
%__________________________________________________________________________   
% PROTOTYPE:
%    plot_flyby(t_dep, t_ga, t_arr, id1, id2, id3)
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
%   figure containing the 3D planet sphere and the powered flyby hyperbolas
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gravitational parameter of the second planet
mu2 = astroConstants(id2 + 10);

% calculate the departure DeltaV and the helicentric arrival velocity at
% the second planet (V_m)
[~, ~, V_m, ~] = single_arc(t_dep, t_ga, id1, id2); 

% calculate the arrival DeltaV and the helicentric departure velocity from
% the second planet (V_p)
[~, ~, ~, V_p] = single_arc(t_ga, t_arr, id2, id3);

% calculate the hyperbolic arcs inside the SOI
[rp, dV_ga, delta_m, delta_p, a_m, a_p, e_m, e_p] =...
    GApow(V_m, V_p, id2, t_ga);

% velocity of second planet at the moment of the flyby
mjd2000_flyby = t_ga; 
[kep2, ksun] = uplanet(mjd2000_flyby, id2);

[R_2, V_2] = kep2car(kep2, ksun);
R_2 = norm(R_2);

% relative velocity of the s/c when approaching the second planet
v_inf_m = V_m - V_2;

% relative velocity of the s/c when leaving the second planet
v_inf_p = V_p - V_2;


% velocity at the pericenter of the first hyperbola


% u vector perpendicular to the orbital plane, used in the findplane
% function
u = cross(v_inf_m, v_inf_p)/norm(cross(v_inf_m, v_inf_p));

% find the plane of the first hyperbola, it will be the same as the second
[i, OM, om] = findplane(u, v_inf_p, delta_m);

% keplerian elements of the two hyperbolic arcs
kep_m = [a_m, e_m, i, OM, om, 0];
kep_p = [a_p, e_p, i, OM, om, 0];

% limits for the two true anomalies of the two hyperbolic arcs
th_inf_m = acos(-1/e_m);
th_inf_p = acos(-1/e_p);

% vectors of true anomalies to be plotted
it = 1e3;
th_m = linspace(-th_inf_m + pi/16, 0, it);
th_p = linspace(0 , th_inf_p - pi/16, it);

% calculate cartesian position of the orbits' cartesian coordinates
r_m = zeros(3,it); r_p = zeros(3,it);
for ii = 1 : it
    kep_m(6) = th_m(ii);
    kep_p(6) = th_p(ii);
    
    [r_m(:, ii), ~] = kep2car(kep_m, mu2);
    [r_p(:, ii), ~] = kep2car(kep_p, mu2);
    
end

% Find 2 points of each asymptote
kep_m_asymptote_inf    = kep_m;
kep_m_asymptote_inf(6) = -th_inf_m + 1e-5;
[r_m_asymptote_inf, ~] = kep2car(kep_m_asymptote_inf, mu2);
r_m_asymptote_0 = r_m(:, end)*e_m;
r_m_asymptote = [r_m_asymptote_0, r_m_asymptote_inf'];

kep_p_asymptote_inf    = kep_p;
kep_p_asymptote_inf(6) = th_inf_p - 1e-5;
[r_p_asymptote_inf, ~] = kep2car(kep_p_asymptote_inf, mu2);
r_p_asymptote_0 = r_p(:, 1)*e_p;
r_p_asymptote = [r_p_asymptote_0, r_p_asymptote_inf'];

% Incoming arc:
plot3(r_m(1,:), r_m(2,:), r_m(3,:), 'b', 'Linewidth', 2)
hold on

% Outgoing arc:
plot3(r_p(1,:), r_p(2,:), r_p(3,:), 'r', 'Linewidth', 2)

% Apse Line
plot3([r_m(1, end)*e_m 0], [r_m(2, end)*e_m 0], [r_m(3, end)*e_m 0], 'k--', 'Linewidth', 1)

axis equal
axis manual

% Incoming asymptote:
plot3(r_m_asymptote(1,:), r_m_asymptote(2,:), r_m_asymptote(3,:), 'g', 'Linewidth', 1)
% Outgoing asymptote:
plot3(r_p_asymptote(1,:), r_p_asymptote(2,:), r_p_asymptote(3,:), 'm', 'Linewidth', 1)

% plot planet 
plot_planet_sphere([],id2);

grid on
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

% plot the DeltaV required by the s/c
% text(r_p(1,1)*1.3, r_p(2,1)*1.3, r_p(3,1)*1.3, ['\DeltaV_{s/c} = ',...
%     num2str(dV_ga), ' km/s'],'color', 'k',...
%     'FontSize',16,'FontName','CMU serif','FontWeight','bold');


leg =legend('Incoming Arc','Outgoing Arc','Apse Line',...
    'Incoming Asymptote','Outgoing Asymptote');
leg.FontName = 'CMU Serif';
leg.FontSize = 18;



ax = gca;

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'X [km]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Y [km]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.String = 'Z [km]';
ax.ZLabel.FontSize = 24;
ax.ZAxis.FontName= 'CMU serif';
ax.ZAxis.FontSize= 16;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % OPTIONAL: plot velocity triangles 
% figure()
% % Planet velocity
% plot3([0 V_2(1)],[0 V_2(2)],[0 V_2(3)],'lineWidth',3) 
% hold on
% % Incoming heliocentric velocity
% plot3([0 V_m(1)],[0 V_m(2)],[0 V_m(3)],'lineWidth',3) 
% % Outgoing heliocentric velocity
% plot3([0 V_p(1)],[0 V_p(2)],[0 V_p(3)],'lineWidth',3) 
% 
% % Incoming planet-centric velocity
% plot3([V_2(1) v_inf_m(1)+V_2(1)],...
%     [V_2(2) v_inf_m(2)+V_2(2)],[V_2(3) v_inf_m(3)+V_2(3)],'lineWidth',3)
% % Outgoing planet-centric velocity
% plot3([V_2(1) v_inf_p(1)+V_2(1)],...
%     [V_2(2) v_inf_p(2)+V_2(2)],[V_2(3) v_inf_p(3)+V_2(3)],'lineWidth',3)
% 
% legend('V_{P}','V^{-}','V^{+}','v^{-}_{inf}','v^{+}_{inf}')
% 
% ax = gca;
% 
% ax.XLabel.String = 'V_X [km/s]';
% ax.XLabel.FontSize = 24;
% ax.XAxis.FontName= 'CMU serif';
% ax.XAxis.FontSize= 16;
% 
% ax.YLabel.String = 'V_Y [km/s]';
% ax.YLabel.FontSize = 24;
% ax.YAxis.FontName= 'CMU serif';
% ax.YAxis.FontSize= 16;
% 
% ax.ZLabel.String = 'V_Z [km/s]';
% ax.ZLabel.FontSize = 24;
% ax.ZAxis.FontName= 'CMU serif';
% ax.ZAxis.FontSize= 16;
% pbaspect([1 2 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: print flyby parameters
% quiver3(0,0,0,v_inf_m(1)*3e3,v_inf_m(2)*3e3,v_inf_m(3)*3e3)
% quiver3(0,0,0,v_inf_p(1)*3e3,v_inf_p(2)*3e3,v_inf_p(3)*3e3)
fprintf('\nFLYBY DATA:\n')

fprintf('\nIncoming hyperbolic arc:')
fprintf('\n a    = %.2f km',kep_m(1))
fprintf('\n e    = %1.4f ',kep_m(2))
fprintf('\n i    = %.2f deg',rad2deg(kep_m(3)))
fprintf('\n OM   = %.2f deg',rad2deg(kep_m(4)))
fprintf('\n om   = %.2f deg',rad2deg(kep_m(5)))
fprintf('\n delta   = %.2f deg',rad2deg(delta_m))
fprintf('\n Delta   = %.2f km',-a_m/tan(delta_m/2))
fprintf('\n\n th_inf = %.2f deg',rad2deg(th_inf_m))

fprintf('\n\nOutgoing hyperbolic arc:')
fprintf('\n a    = %.2f km',kep_p(1))
fprintf('\n e    = %1.4f ',kep_p(2))
fprintf('\n i    = %.2f deg',rad2deg(kep_p(3)))
fprintf('\n OM   = %.2f deg',rad2deg(kep_p(4)))
fprintf('\n om   = %.2f deg',rad2deg(kep_p(5)))
fprintf('\n delta   = %.2f deg',rad2deg(delta_p))
fprintf('\n Delta   = %.2f km',-a_p/tan(delta_p/2))
fprintf('\n\n th_inf = %.2f deg',rad2deg(th_inf_p))

fprintf('\n\n Pericentre Altitude = %.2f km',rp - astroConstants(20+id2))
fprintf('\n Total DeltaV = %.4f km/s',norm(v_inf_p-v_inf_m))
fprintf('\n Velocity at pericentre = %.2f km/s',sqrt(mu2*(2/rp-1/a_m)))
fprintf('\n (taken before powered DeltaV)')


% find time spent inside SOI
r_SOI = R_2 * (astroConstants(10+id2)/ksun)^(2/5);
th_SOI_m = 2*pi - acos(1/e_m * (a_m*(1-e_m^2) / r_SOI - 1));
th_SOI_p = acos(1/e_p * (a_p*(1-e_p^2) / r_SOI - 1));
F_SOI_m = 2*atanh(sqrt((e_m-1)/(e_m+1))*tan(th_SOI_m/2));
F_SOI_p = 2*atanh(sqrt((e_p-1)/(e_p+1))*tan(th_SOI_p/2));
Dt_m = abs(sqrt(-a_m^3 / mu2) * (e_m*sinh(F_SOI_m) - F_SOI_m));
Dt_p = abs(sqrt(-a_p^3 / mu2) * (e_p*sinh(F_SOI_p) - F_SOI_p));
Dt = (Dt_m + Dt_p)/3600;
fprintf('\n Time spent in SOI = %.2f h\n\n',Dt)


