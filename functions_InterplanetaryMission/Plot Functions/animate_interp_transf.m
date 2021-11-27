function animate_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3,s)
% ANIMATE INTERPLANETARY TRANSFER: plots a 3D animated graph of the two 
%   interplanetary legs of the three-planet transfer, together with the 
%   three planet's orbits. 
%   The current date will be plotted on the graph.
%   At the bottom of the function there's a part that may be uncommented by
%   the user in order to have extra data printed in the command window
%__________________________________________________________________________   
% PROTOTYPE:
%    animate_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3,s)
% 
% INPUT:
%   t_dep[1]        departure time in MJD2000                       [days]
%   dt_ga[1]        gravity assist time in MJD2000                  [days]
%   dt_arr[1]       arrival time in MJD2000                         [days]
%   id1[1]          departure planet identifier                     [-]
%   id2[1]          gravity assist planet identifier                [-]
%   id3[1]          arrival planet identifier                       [-]
%   s[1]            speed of the animation                          [-]
%
% OUTPUT:
%   figure containing the 3D animated plot
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The plot is in astronomical units:
AU = astroConstants(2);

% Find useful initial parameters
[kep1,muSun] = uplanet(t_dep, id1);
[r0_1,v1] = kep2car(kep1, muSun);
[kep2,muSun] = uplanet(t_dep, id2);
[r0_2,v2] = kep2car(kep2, muSun);
[kep3,muSun] = uplanet(t_dep, id3);
[r0_3,v3] = kep2car(kep3, muSun);

[kep2,muSun] = uplanet(t_ga , id2);
[r2_ga ,~] = kep2car(kep2, muSun);
[kep3,muSun] = uplanet(t_arr, id3);
[r3_arr,~] = kep2car(kep3, muSun);

Dt1 = (t_ga - t_dep)*24*3600;
[~, ~, ~, ~, V_p1, ~, ~, ~] = lambertMR( r0_1, r2_ga, Dt1, muSun, 0, 0, 0, 2 );
Dt2 = (t_arr - t_ga)*24*3600;
[~, ~, ~, ~, V_p2, ~, ~, ~] = lambertMR( r2_ga, r3_arr, Dt2, muSun, 0, 0, 0, 2 );

% define integration parameters
tspan = linspace(t_dep*24*3600,t_arr*24*3600,1e5);
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
y0_1  = [r0_1'  ; v1'  ];
y0_2  = [r0_2'  ; v2'  ];
y0_3  = [r0_3'  ; v3'  ];
y0_T1 = [r0_1'  ; V_p1'];
y0_T2 = [r2_ga' ; V_p2'];

% integrate planet and s/c orbits
[~,y_1] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0_1,options);
r_1 = y_1(:,1:3)/AU;

[~,y_2] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0_2,options);
r_2 = y_2(:,1:3)/AU;

[~,y_3] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0_3,options);
r_3 = y_3(:,1:3)/AU;

[~,y_T1] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0_T1,options);
r_T1 = y_T1(:,1:3)/AU;

[~,y_T2] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
               tspan,y0_T2,options);
r_T2 = y_T2(:,1:3)/AU;


% Plot planet obits and setup figure
plot_planet_orbit_k(t_dep,t_dep,id1)
hold on
set(gca,'Color','k')
set(gcf,'color','k')
axis equal
grid on
plot_planet_orbit_k(t_dep,t_dep,id2)
plot_planet_orbit_k(t_dep,t_dep,id3)
xlim([-1.05 1.05])
ylim([-1.05 1.05])
zlim([-0.1 0.1])
fig = gcf;
fig.Position = 1.0e+03 *[0.3497    0.3070    1.8833    1.0293];

% Plot the Sun
plot_planet_sphere(t_dep,0,6e6)

% Name axes
ax = gca;


ax.GridColor = 'w';
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
% Initialize animated elements
orbit_1    = animatedline(r_1(1,1),r_1(1,2),r_1(1,3),'color','w','linew',2);
orbit_2    = animatedline(r_2(1,1),r_2(1,2),r_2(1,3),'color','w','linew',2);
orbit_3    = animatedline(r_3(1,1),r_3(1,2),r_3(1,3),'color','w','linew',2);

planet_1 = plot_planet_sphere_handle([],id1,3e-2);
XData_1=planet_1.XData; YData_1=planet_1.YData; ZData_1=planet_1.ZData;
planet_2 = plot_planet_sphere_handle([],id2,3e-2);
XData_2=planet_2.XData; YData_2=planet_2.YData; ZData_2=planet_2.ZData;
planet_2_fb = plot_planet_sphere_handle([],id2,3e-2);
XData_2_fb=planet_2_fb.XData; YData_2_fb=planet_2_fb.YData; ZData_2_fb=planet_2_fb.ZData;
planet_3 = plot_planet_sphere_handle([],id3,3e-2);
XData_3=planet_3.XData; YData_3=planet_3.YData; ZData_3=planet_3.ZData;

transfer_1 = animatedline(r_T1(1,1),r_T1(1,2),r_T1(1,3),'color','g','linew',3);
transfer_2 = animatedline(r_T2(1,1),r_T2(1,2),r_T2(1,3),'color','r','linew',3);

date_vect = mjd20002date(t_dep);
date = text(1,1,0,[num2str(date_vect(1)) '/' num2str(date_vect(2))...
     '/' num2str(date_vect(3))]);
 
date.FontName = 'CMU Serif';
date.FontSize = 20;
date.Color = 'w';

sc = plot3(r_T1(1,1),r_T1(1,2),r_T1(1,3),'r.','MarkerSize',20);
k = 0;
flag = false;
ax.View = [-74.6575    9.8207];
% Animate using a for cycle
for ii = 1 : s : length(tspan)
    % check if flyby has already happened
    if tspan(ii) > t_ga*24*3600
        afterflyby = true;
    else
        afterflyby = false;
    end
    
    % update planet orbit
    addpoints(orbit_1,r_1(ii,1),r_1(ii,2),r_1(ii,3))
    addpoints(orbit_2,r_2(ii,1),r_2(ii,2),r_2(ii,3))
    addpoints(orbit_3,r_3(ii,1),r_3(ii,2),r_3(ii,3))
    
    % update s/c orbit and position
    if afterflyby == false
        addpoints(transfer_1,r_T1(ii,1),r_T1(ii,2),r_T1(ii,3))
        sc.XData=r_T1(ii,1); sc.YData=r_T1(ii,2); sc.ZData=r_T1(ii,3);
        

        k = k+1;
    elseif afterflyby == true
        addpoints(transfer_2,r_T2(ii-k*s,1),r_T2(ii-k*s,2),r_T2(ii-k*s,3))
        sc.XData=r_T2(ii-k*s,1); sc.YData=r_T2(ii-k*s,2); sc.ZData=r_T2(ii-k*s,3);
        

        if flag == false
            % Leave venus at flyby point
            planet_2_fb.XData = XData_2_fb + r_2(ii,1);
            planet_2_fb.YData = YData_2_fb + r_2(ii,2);
            planet_2_fb.ZData = ZData_2_fb + r_2(ii,3);
            flag = true;
        end
    end
    
    % update planet sphere position
    planet_1.XData = XData_1 + r_1(ii,1);
    planet_1.YData = YData_1 + r_1(ii,2);
    planet_1.ZData = ZData_1 + r_1(ii,3);

    planet_2.XData = XData_2 + r_2(ii,1);
    planet_2.YData = YData_2 + r_2(ii,2);
    planet_2.ZData = ZData_2 + r_2(ii,3);

    planet_3.XData = XData_3 + r_3(ii,1);
    planet_3.YData = YData_3 + r_3(ii,2);
    planet_3.ZData = ZData_3 + r_3(ii,3);
    
    % update date
    mjd2000 = t_dep + (tspan(ii)-tspan(1)) / (3600*24);
    date_vect = mjd20002date(mjd2000);
    date.String = [num2str(date_vect(1)) '/' num2str(date_vect(2))...
     '/' num2str(date_vect(3))];
    
    drawnow
    pause(eps)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




