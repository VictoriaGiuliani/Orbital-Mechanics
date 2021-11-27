function plot_neighbourhood(t_dep,t_ga_n,t_arr_n,id1,id2,id3,days_ga,days_arr,it)
% PLOT NEIGHBOURHOOD: plots three contour graphs around certain time 
%   instants specified by the user. Time of departure is kept constant,
%   while time of flights for gravity assist and arrival are changed to
%   form the two axes of the contour plots.
%   The Z axis values are the DeltaVs associated to the transfers with the
%   aformentioned time of departure, gravity assist and arrival.
%   White zones on the plots mean that the transfer is not physically
%   feasible (e.g. satellite crashes on flyby planet).
%__________________________________________________________________________   
% PROTOTYPES:
%    plot_neighbourhood(t_dep_n,dt_ga_n,dt_arr_n,id1,id2,id3,days_ga,days_arr,it)
% or plot_neighbourhood(t_dep_n,dt_ga_n,dt_arr_n,id1,id2,id3,days_ga,days_arr)
% or plot_neighbourhood(t_dep_n,dt_ga_n,dt_arr_n,id1,id2,id3)
% 
% INPUT:
%   t_dep[1]        departure time in MJD2000                       [days]
%   t_ga_n[1]       gravity assist time in MJD2000, around which the
%                   neighbourhood will be found                     [days]
%   t_ga_n[1]       arrival time in MJD2000, around which the
%                   neighbourhood will be found                     [days]
%   id1[1]          departure planet identifier                     [-]
%   id2[1]          gravity assist planet identifier                [-]
%   id3[1]          arrival planet identifier                       [-]
%   days_ga[1]      xlim of the contour plot, specifies how large
%                   the neighbourhood is in the dt_ga direction     [days]
%   days_ga[1]      ylim of the contour plot, specifies how large
%                   the neighbourhood is in the dt_arr direction    [days]
%   it[1]           number of iterations performed both in the
%                   dt_ga timespan and in the dt_arr timespan, this
%                   means that it^2 transfers are calculated        [-]
%
% OUTPUT:
%   figure containing contour plots of the total DeltaV, the DeltaV
%       spent in the gravity assist and the sum of departure and
%       arrival DeltaVs
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check number of inputs, define default values not set by the user
if nargin == 6
    days_ga  = 50;
    days_arr = 50;
    it = 100;
elseif nargin == 8
    it = 100;
elseif nargin ~= 9
    error('Wrong number of inputs')
end

% Find transfer times
dt_ga_n  = t_ga_n  - t_dep;
dt_arr_n = t_arr_n - t_ga_n;

% Define arrays for the for cycle

dt_ga  = linspace(max([20 dt_ga_n - days_ga]), dt_ga_n + days_ga,it);
dt_arr = linspace(max([20 dt_arr_n - days_arr]), dt_arr_n + days_arr,it);
dV_tot = zeros(it,it);
dV_dep = zeros(it,it);
dV_arr = zeros(it,it);
dV_ga = zeros(it,it);

% Perform all it^2 transfers, t_dep is kept constant while t_ga and t_arr
% change
for j = 1 : length(dt_ga)
    for k = 1 : length(dt_arr)

        [dV_tot(j, k),dV_dep(j, k),dV_arr(j, k),dV_ga(j, k)] =...
            GA_interp_transf(t_dep, t_dep+dt_ga(j),...
            t_dep+dt_ga(j)+dt_arr(k), id1, id2, id3);
    end
end

% Perform contour plots
%__________________________________________________________________________
subplot(2,3,[1 2 4 5])
contourf(dt_ga,dt_arr,dV_tot')
hold on
global mindV
if ~isempty(mindV)
    mindV = [mindV min(min(dV_tot))];
end
[i,j] = find(dV_tot == min(min(dV_tot)));
scatter(dt_ga(i),dt_arr(j),'or','filled')

ax = gca;

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$\Delta t_{ga}$ [days]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\Delta t_{arr}$ [days]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

ax.Title.FontSize = 24;
ax.Title.FontName = 'CMU Serif';
date_dep = datestr(mjd20002date(t_dep),1);
ax.Title.String = sprintf('Total \\DeltaV, departure on %s',date_dep);

cb = colorbar;
cb.FontName = 'CMU Serif';
cb.FontSize = 16;
cb.Label.String = '\DeltaV [km/s]';
%__________________________________________________________________________
subplot(2,3,3)
contourf(dt_ga,dt_arr,dV_ga')
hold on
scatter(dt_ga(i),dt_arr(j),'or','filled')
ax = gca;

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$\Delta t_{ga}$ [days]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\Delta t_{arr}$ [days]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

ax.Title.FontSize = 24;
ax.Title.FontName = 'CMU Serif';
ax.Title.String = '\DeltaV_{ga}';

cb = colorbar;
cb.FontName = 'CMU Serif';
cb.FontSize = 16;
cb.Label.String = '\DeltaV [km/s]';
%__________________________________________________________________________
subplot(2,3,6)
contourf(dt_ga,dt_arr,dV_dep'+dV_arr')
hold on
scatter(dt_ga(i),dt_arr(j),'or','filled')
ax = gca;

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$\Delta t_{ga}$ [days]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\Delta t_{arr}$ [days]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

ax.Title.FontSize = 24;
ax.Title.FontName = 'CMU Serif';
ax.Title.String = '\DeltaV_{dep} + \DeltaV_{arr}';

cb = colorbar;
cb.FontName = 'CMU Serif';
cb.FontSize = 16;
cb.Label.String = '\DeltaV [km/s]';
%__________________________________________________________________________
% Optional: 3D plots of the contour plots
% figure(3)
% [X,Y]=meshgrid(dt_ga,dt_arr);
% surf(X,Y,dV_tot')
% figure(4)
% surf(X,Y,dV_ga')
% figure(5)
% surf(X,Y,dV_dep'+dV_arr')