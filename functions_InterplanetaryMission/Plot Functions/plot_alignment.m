function plot_alignment(t_dep_n,dt_ga_n,dt_arr_n,id1,id2,id3)
% PLOT OF ALIGNMENT: plots the transfer legs around the Mercury-Venus
%   alignment. (Used only for this specific case).
%__________________________________________________________________________   
% PROTOTYPE:
%    plot_alignment(t_dep_n,dt_ga_n,dt_arr_n,id1,id2,id3)
%
% INPUT:
%   t_dep[1]  time in MJD2000, for the departure of the s/c.        [days]
%             (at "alignment")
%   dt_ga[1]  days between departure and flyby                      [days]
%             (at "alignment")
%   dt_arr[1] days between flyby and arrival                        [days]
%             (at "alignment")
%   id1[1]    departure planet identifier                           [-]
%   id2[1]    gravity assist planet identifier                      [-]
%   id3[1]    arrival planet identifier                             [-]
%
% OUTPUT:
%   figure containing the various transfers near the alignment one.
%   figure containing the departure DeltaV wrt days from alignment.
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% useful parameters:
days = 20; % half width of neighbourhood to be plotted (centre at alignment)
it = 30; % number of transfers to plot

% set time (flyby is a vector)
t_dep = t_dep_n;
t_ga  = t_dep + linspace( dt_ga_n-days , dt_ga_n+days , it*10);
t_arr = t_dep + dt_ga_n + dt_arr_n;

dV_tot = zeros(it*10,1); dV_dep = zeros(it*10,1);
dV_arr = zeros(it*10,1); dV_ga  = zeros(it*10,1);

% find all DeltaVs
for i = 1:it*10   
    [dV_tot(i),dV_dep(i),dV_arr(i),dV_ga(i)] =...
              GA_interp_transf(t_dep, t_ga(i), t_arr, id1, id2, id3);    
end



% The plot is in astronomical units:
AU = astroConstants(2);

% calculate transfers in for loop
for i = 1:it
    % Calculate the first leg of the transfer using the lambert solver and
    % integrating it in cartesian form
    [kep1,muSun] = uplanet(t_dep, id1);
    [r1,~] = kep2car(kep1, muSun);

    [kep2,muSun] = uplanet(t_ga(i*10), id2);
    [r2,~] = kep2car(kep2, muSun);
    Dt = (t_ga(i*10) - t_dep)*24*3600;
    [~, ~, ~, ~, V_p, ~, ~, ~] = lambertMR( r1, r2, Dt, muSun, 0, 0, 0, 2 );

    tspan = linspace(t_dep*24*3600,t_ga(i*10)*24*3600,1000);
    y0 = [r1' ; V_p'];

    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
    [~,y1] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
                   tspan,y0,options);

               y1 = y1/AU;
    % Plot first leg           
    plot3(y1(:,1),y1(:,2),y1(:,3),'Color',[i/it 0.5 1-i/it],'LineWidth',2)
    hold on


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the second leg of the transfer using the lambert solver and
    % integrating it in cartesian form
    [kep1,muSun] = uplanet(t_ga(i*10), id2);
    [r1,~] = kep2car(kep1, muSun);

    [kep2,muSun] = uplanet(t_arr, id3);
    [r2,~] = kep2car(kep2, muSun);
    Dt = (t_arr - t_ga(i*10))*24*3600;
    [~, ~, ~, ~, V_p, ~, ~, ~] = lambertMR( r1, r2, Dt, muSun, 0, 0, 0, 2 );

    tspan = linspace(t_ga(i*10)*24*3600,t_arr*24*3600,1000);
    y0 = [r1' ; V_p'];

    options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
    [~,y2] = ode113(@(t,y) ODE_two_body(t,y,muSun),...
                   tspan,y0,options);

               y2 = y2/AU;
    % Plot second leg           
    plot3(y2(:,1),y2(:,2),y2(:,3),'b','LineWidth',1)
    
    
    plot_planet_orbit(t_dep,t_ga(i*10) ,id2)
    
end
cb = colorbar;

cb.Ticks = [-days:4:days]';
for ii = 1:length(cb.Ticks)
    cb.TickLabels{ii} = num2str(cb.Ticks(ii));
end
% attention, colormap must be personalized by user
cb.Limits = [-days days];
cb.FontName = 'CMU Serif';
cb.FontSize = 16;
cb.Label.String = 'Days from Alignment';
% Plot the Sun
plot_planet_sphere(t_dep,0,6e6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the orbits of the three planets, highliting the portions spent
% between departure and closest approach to s/c
plot_planet_orbit(t_dep,t_dep,id1)
plot_planet_sphere(t_dep,id1,3e6)

plot_planet_sphere(t_dep,id2,3e6)


plot_planet_orbit(t_dep,t_arr,id3)
plot_planet_sphere(t_dep,id3,3e6)
plot_planet_sphere(t_arr,id3,3e6)

col = 'w';
axis equal; pbaspect([1 1 1]); set(gca,'Color',col);
grid on


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
figure(2)

% figure containing the departure DeltaV wrt days from alignment.
plot(linspace(-days,days,it*10),dV_dep,'LineWidth',4);
ax = gca;
grid on
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Days from Alignment';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.String = '\DeltaV_{dep} [km/s]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;
end

