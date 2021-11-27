function plot_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3)
% PLOT INTERPLANETARY TRANSFER: plots a 3D graph of the two interplanetary
%   legs of the three-planet transfer, together with the three planet's
%   orbits. The orbits of the second and third planets will have highlited
%   parts, representing the planets' flown trajectory from the time of 
%   departure to the time of closest approach to the s/c.
%   DeltaVs will be printed near each manoeuvre point, the total DeltaV is
%   printed near the Sun.
%   At the bottom of the function there's a part that may be uncommented by
%   the user in order to have extra data printed in the command window
%__________________________________________________________________________   
% PROTOTYPE:
%    plot_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3)
% 
% INPUT:
%   t_dep[1]        departure time in MJD2000                       [days]
%   dt_ga[1]        gravity assist time in MJD2000                  [days]
%   dt_arr[1]       arrival time in MJD2000                         [days]
%   id1[1]          departure planet identifier                     [-]
%   id2[1]          gravity assist planet identifier                [-]
%   id3[1]          arrival planet identifier                       [-]
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

% Calculate the first leg of the transfer using the lambert solver and
% integrating it in cartesian form
[kep1,muSun] = uplanet(t_dep, id1);
[r1,v1] = kep2car(kep1, muSun);

[kep2,muSun] = uplanet(t_ga, id2);
[r2,~] = kep2car(kep2, muSun);
Dt1 = (t_ga - t_dep)*24*3600;
[~, ~, ~, ~, V_p, ~, ~, Dth1] = lambertMR( r1, r2, Dt1, muSun, 0, 0, 0, 2 );

% Calculate keplerian parameters of first leg
kep1_leg = car2kep (r1, V_p, muSun);
th_end_1 = wrapTo2Pi(kep1_leg(6) + Dth1);
% Calculate dihedral angle between planet and transfer's plane
h_planet1 = cross(r1,v1);
h_transfer1 = cross(r1,V_p);
dihedral_1 = acos(dot(h_planet1,h_transfer1)/...
    (norm(h_planet1)*norm(h_transfer1)));


tspan = linspace(t_dep*24*3600,t_ga*24*3600,1000);
y0 = [r1' ; V_p'];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[~,y1] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0,options);

           y1 = y1/AU;
% Plot first leg           
plot3(y1(:,1),y1(:,2),y1(:,3),'LineWidth',4)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the second leg of the transfer using the lambert solver and
% integrating it in cartesian form
[kep1,muSun] = uplanet(t_ga, id2);
[r1,~] = kep2car(kep1, muSun);

[kep2,muSun] = uplanet(t_arr, id3);
[r2,v2] = kep2car(kep2, muSun);
Dt2 = (t_arr - t_ga)*24*3600;
[~, ~, ~, ~, V_p, ~, ~, Dth2] = lambertMR( r1, r2, Dt2, muSun, 0, 0, 0, 2 );

% Calculate keplerian parameters of second leg
kep2_leg = car2kep (r1, V_p, muSun);
th_end_2 = wrapTo2Pi(kep2_leg(6) + Dth2);
% Calculate dihedral angle between first and second leg planes
h_transfer2 = cross(r1,V_p);
dihedral_2 = acos(dot(h_transfer1,h_transfer2)/...
    (norm(h_transfer1)*norm(h_transfer2)));
% Calculate dihedral angle between second leg and third planet's planes
h_planet3 = cross(r2,v2);
dihedral_3 = acos(dot(h_planet3,h_transfer2)/...
    (norm(h_planet3)*norm(h_transfer2)));


tspan = linspace(t_ga*24*3600,t_arr*24*3600,1000);
y0 = [r1' ; V_p'];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[~,y2] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0,options);
           
           y2 = y2/AU;
% Plot second leg           
plot3(y2(:,1),y2(:,2),y2(:,3),'LineWidth',4)

% Plot the Sun
plot_planet_sphere(t_dep,0,6e6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the orbits of the three planets, highliting the portions spent
% between departure and closest approach to s/c
plot_planet_orbit(t_dep,t_dep,id1)
plot_planet_sphere(t_dep,id1,3e6)

plot_planet_orbit(t_dep,t_ga ,id2)
plot_planet_sphere(t_dep,id2,3e6)
plot_planet_sphere(t_ga ,id2,3e6)

plot_planet_orbit(t_dep,t_arr,id3)
plot_planet_sphere(t_dep,id3,3e6)
plot_planet_sphere(t_arr,id3,3e6)

col = 'w';
if col == 'w'
    col_opp = 'k';
else
    col_opp = 'w';
end
axis equal; pbaspect([1 1 1]); set(gca,'Color',col);
grid on

% Calculate and display the DeltaVs at each manoeuvre point, the total
% DeltaV is displayed near the sun
[dV_tot,dV_dep,dV_arr,dV_ga] =...
                GA_interp_transf(t_dep,t_ga,t_arr, id1, id2, id3);
            
text(0, -0.1, 0, ['\DeltaV_{tot} = ',...
    num2str(dV_tot), ' km/s'],'color', col_opp,...
    'FontSize',16,'FontName','CMU serif','FontWeight','bold');
text(y1(1,1), y1(1,2)-0.1, y1(1,3), ['\DeltaV_{dep} = ',...
    num2str(dV_dep), ' km/s'],'color', col_opp,...
    'FontSize',16,'FontName','CMU serif','FontWeight','bold');
text(y1(end,1), y1(end,2)-0.1, y1(end,3), ['\DeltaV_{ga} = ',...
    num2str(dV_ga), ' km/s'],'color', col_opp,...
    'FontSize',16,'FontName','CMU serif','FontWeight','bold');
text(y2(end,1), y2(end,2)-0.1, y2(end,3), ['\DeltaV_{arr} = ',...
    num2str(dV_arr), ' km/s'],'color', col_opp,...
    'FontSize',16,'FontName','CMU serif','FontWeight','bold');


% Name axes
ax = gca;

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'X [AU]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Y [AU]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.String = 'Z [AU]';
ax.ZLabel.FontSize = 24;
ax.ZAxis.FontName= 'CMU serif';
ax.ZAxis.FontSize= 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL: print leg parameters on command window
fprintf('\nFirst leg Keplerian parameters:')
fprintf('\n a    = %.2f km',kep1_leg(1))
fprintf('\n e    = %1.4f ',kep1_leg(2))
fprintf('\n i    = %.2f deg',rad2deg(kep1_leg(3)))
fprintf('\n OM   = %.2f deg',rad2deg(kep1_leg(4)))
fprintf('\n om   = %.2f deg',rad2deg(kep1_leg(5)))
fprintf('\n th_0 = %.2f deg',rad2deg(kep1_leg(6)))
fprintf('\n th_f = %.2f deg\n',rad2deg(th_end_1))
fprintf('\n Dt   = %.2f days\n',Dt1/24/3600)

fprintf('\nSecond leg Keplerian parameters:')
fprintf('\n a    = %.2f km',kep2_leg(1))
fprintf('\n e    = %1.4f ',kep2_leg(2))
fprintf('\n i    = %.2f deg',rad2deg(kep2_leg(3)))
fprintf('\n OM   = %.2f deg',rad2deg(kep2_leg(4)))
fprintf('\n om   = %.2f deg',rad2deg(kep2_leg(5)))
fprintf('\n th_0 = %.2f deg',rad2deg(kep2_leg(6)))
fprintf('\n th_f = %.2f deg\n',rad2deg(th_end_2))
fprintf('\n Dt   = %.2f days\n',Dt2/24/3600)

fprintf('\nDihedral angles between planes\n at each manoeuvre:')
fprintf('\n delta M-1 = %.2f deg',rad2deg(dihedral_1))
fprintf('\n delta 1-2 = %.2f deg',rad2deg(dihedral_2))
fprintf('\n delta 2-E = %.2f deg\n',rad2deg(dihedral_3))
