%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INTERPLANETARY EXPLORER SCRIPT
%
%   For faster use, execute just the first part of this code, then load 
%   one of the saved workspaces in this folder.
%   Example: 40_yrs_2080_it.mat means 40 years departure window with 2080
%            iterations on it.
%   This allows the user to skip the computationally expensive FOR SOLUTION
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% Adding to path necessary functions
addpath('Basic Functions','Planet Textures','Plot Functions')

% Input the planets (in this case Mercury, Venus, Earth)
id1 = 1;
id2 = 2;
id3 = 3;
% Height of the second planet's atmosphere (used in flyby)
global h_atm 
h_atm = 150;


% Find the lower and upper bounds of the search region in the for cycle
[dt_ga_lower ,dt_ga_upper ] = twindow(id1,id2,0.5,2);
[dt_arr_lower,dt_arr_upper] = twindow(id2,id3,0.5,2);

% Departure window and conversion to MJD2000, for the upper boundary of the
% departure window dt_ga_lower and dt_arr_lower are subtracted because any
% later transfer wouldn't arrive in time (arrival time window assumed equal
% to the departure time window)
t_dep_win = [2027 12 1 0 0 0; 2067 12 1 0 0 0];
t1_start_MJD2000 = date2mjd2000(t_dep_win(1,1:6));
t1_end_MJD2000   = date2mjd2000(t_dep_win(2,1:6)) - dt_ga_lower - dt_arr_lower;

% Iterations in the for cycle
it_dep = 400;    % for departure window
it_ga  = 20;     % for gravity-assist window
it_arr = 20;     % for arrival window

% Time of departure vector used in the for cycle
t_dep = linspace(t1_start_MJD2000, t1_end_MJD2000, it_dep);

% Initialize arrays for the for cycle
dV_tot = zeros(it_dep,it_ga,it_arr);dV_dep = zeros(it_dep,it_ga,it_arr);
dV_arr = zeros(it_dep,it_ga,it_arr);dV_ga  = zeros(it_dep,it_ga,it_arr);
% Initialize feasible solutions counter
n = 0;



%% FOR SOLUTION
% WARNING: to execute this part in a more stable way, change the conditions
%   on the pericentre radius in the GApow function (put && instead of ||)
%
% Rough grid search of the departure, gravity assist and arrival time
% windows
for i = 1 : length(t_dep)
    % Definition of time of flyby vector in MJD2000
    t_ga = linspace(t_dep(i) + dt_ga_lower, t_dep(i) + dt_ga_upper, it_ga);
    for j = 1 : length(t_ga)
        % Definition of time of arrival vector in MJD2000
        t_arr = linspace(t_ga(j) + dt_arr_lower, t_ga(j) + dt_arr_upper, it_arr);
        for k = 1 : length(t_arr)
            % Calculate interplanetary transfer
            [dV_tot(i, j, k),dV_dep(i, j, k),dV_arr(i, j, k),dV_ga(i, j, k)] =...
                GA_interp_transf(t_dep(i), t_ga(j), t_arr(k), id1, id2, id3);
            
            % If the transfer is feasible, the dep, ga, arr times are
            % recorded
            if ~isnan(dV_tot(i, j, k))
               n = n + 1;
               T(n,:) = [t_dep(i), t_ga(j), t_arr(k)];
            end
        end
    end
end


% initialize vectors for the for cycle
dV_tot1 = zeros(n,1); dV_dep1 = zeros(n,1);
dV_arr1 = zeros(n,1); dV_ga1  = zeros(n,1);
X=[];
% options for the fminunc function
options = optimoptions('fminunc','Display','off');

% initialize computation time measurement
tic
for i = 1 : n
    % call the optimization function with the already calculated feasible 
    % times as initial guess
    x = fminunc(@(x) GA_interp_transf (x(1),x(2),x(3),id1,id2,id3),T(i,:),options);
    % store local minima times
    X = [X;x];
    % store local minima DeltaVs
    [dV_tot1(i),dV_dep1(i),dV_arr1(i),dV_ga1(i)] =...
                GA_interp_transf(x(1),x(2),x(3), id1, id2, id3);
    fprintf('Finding local minima: %d/%d\n',i,n)
end
% save computation time measurement
computationtime = toc;

% save new T matrix with departure time in MJD2000 in the first column,
% time of flight in days until the flyby in the second column, and time of
% flight in days form flyby to arrival in the third column
T(:,1) = X(:,1);
T(:,2) = X(:,2) - X(:,1);
T(:,3) = X(:,3) - X(:,2);

% save workspace
save('data_opt')

%% Local minimums plots
addpath('Basic Functions','Planet Textures','Plot Functions')
close all
max_dV_plot = 30;
T_scatter = T(dV_tot1<max_dV_plot,:);
dV_scatter = dV_tot1(dV_tot1<max_dV_plot);

syn_T = 1.59 * astroConstants(32); % Synodic period in days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first 3D scatter plot
dep = datetime(mjd20002jd(T_scatter(:,1)),'ConvertFrom','juliandate');
scatter3 (dep,T_scatter(:,2),T_scatter(:,3),[],dV_scatter,'filled')
xlim([datetime(mjd20002jd(t1_start_MJD2000),'ConvertFrom','juliandate')...
    datetime(mjd20002jd(t1_end_MJD2000),'ConvertFrom','juliandate')])

pbaspect([3 2 1])

ax = gca;

ax.XLabel.String = 'Departure Date [dd/mm/yyyy]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;
xtickformat('dd/MM/yyyy')


ax.YLabel.String = 'First Transfer Leg ToF [days]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

ax.ZLabel.String = 'Second Transfer Leg ToF [days]';
ax.ZLabel.FontSize = 24;
ax.ZAxis.FontName= 'CMU serif';
ax.ZAxis.FontSize= 16;

cb = colorbar;
cb.FontName = 'CMU Serif';
cb.FontSize = 16;
cb.Label.String = '\DeltaV_{tot} [km/s]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% second 2D scatter plot
figure(2)
scatter (dep,dV_scatter,[],T_scatter(:,2)+T_scatter(:,3),'filled')

xlim([datetime(mjd20002jd(t1_start_MJD2000),'ConvertFrom','juliandate')...
    datetime(mjd20002jd(t1_start_MJD2000 + 10*365.25),'ConvertFrom','juliandate')])
ylim([0 max_dV_plot])

ax = gca;

ax.XLabel.String = 'Departure Date [dd/mm/yyyy]';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;
xtickformat('dd/MM/yyyy')


ax.YLabel.String = '\DeltaV_{tot} [km/s]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

cb = colorbar;
cb.FontName = 'CMU Serif';
cb.FontSize = 16;
cb.Label.String = 'Total Transfer Time [days]';
hold on
for i = 1:floor((t1_end_MJD2000 - t1_start_MJD2000)/syn_T)
    plot([datetime(mjd20002jd(t1_start_MJD2000+syn_T*i),'ConvertFrom','juliandate') ...
        datetime(mjd20002jd(t1_start_MJD2000+syn_T*i),'ConvertFrom','juliandate')]...
        ,[0 max_dV_plot],'k--','lineWidth',1);
end


%% Neighbourhood plots
close all

% data for crescent moon pattern
t_dep_n = 12395;
dt_ga_n  = 136;
dt_arr_n = 163;


plot_neighbourhood(t_dep_n, t_dep_n+dt_ga_n, t_dep_n+dt_ga_n+dt_arr_n,...
    id1,id2,id3,40,50,100)

figure(2)
plot_interp_transf(t_dep_n, dt_ga_n+t_dep_n, dt_ga_n+t_dep_n+dt_arr_n,id1,id2,id3)

%% Alignment plot
close all
% times for Venus and Mercury alignment
t_dep_n = 12395;
dt_ga_n  = 106;
dt_arr_n = 163;

% ecliptic coordinates of the planets
    AU = astroConstants(2);
    [kep,mu_S]= uplanet(t_dep_n + dt_ga_n, id2);
    [r0,~] = kep2car( kep, mu_S );
    fprintf('Venus position at flyby:  %f %f %f\n',r0(1)/AU,r0(2)/AU,r0(3)/AU)
    [kep,mu_S]= uplanet(t_dep_n, id1);
    [r0,~] = kep2car( kep, mu_S );
    fprintf('Mercury initial position: %f %f %f\n',r0(1)/AU,r0(2)/AU,r0(3)/AU)

% plot alignment
plot_alignment(t_dep_n,dt_ga_n,dt_arr_n,id1,id2,id3)


%% GLOBAL OPTIMUM PLOTS
close all
clc
i = ind2sub(size(dV_tot1),find(dV_tot1 == min(min(min(dV_tot1)))));

t_dep = T(i,1);
t_ga  = T(i,2) + t_dep;
t_arr = T(i,3) + t_ga ;

plot_neighbourhood(t_dep,t_ga,t_arr,id1,id2,id3,50,50,100)
figure
plot_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3)
figure
plot_flyby(t_dep,t_ga,t_arr, id1, id2, id3)
% figure
% animate_interp_transf(t_dep,t_ga,t_arr,id1,id2,id3,100)

date_dep = mjd20002date(t_dep);
date_ga  = mjd20002date(t_ga);
date_arr = mjd20002date(t_arr);
fprintf('\nDeparture Date: %d / %d / %d\n',date_dep(1),date_dep(2),date_dep(3));
fprintf('Flyby Date:     %d / %d / %d\n',date_ga(1),date_ga(2),date_ga(3));
fprintf('Arrival Date:   %d / %d / %d\n',date_arr(1),date_arr(2),date_arr(3));

fprintf('Total Transfer time: %3.0f days\n',t_arr-t_dep);



%% DEPARTURE DATE DEPENDENCE

% quite inefficient, recommend not running, results saved in current folder
% as dep_neig.mat
global mindV 
mindV = 0;
dep_search = -7:0.5:7; % days
for j = dep_search
    t_dep = T(i,1) + j;
    t_ga  = T(i,2) + t_dep;
    t_arr = T(i,3) + t_ga ;
    plot_neighbourhood(t_dep,t_ga,t_arr,id1,id2,id3,50,50,100)
    close all
    j
end
mindV(1)=[];
plot(dep_search,mindV,'lineWidth',3)
grid on
ax = gca;

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Days from optimal departure date';
ax.XLabel.FontSize = 24;
ax.XAxis.FontName= 'CMU serif';
ax.XAxis.FontSize= 16;

ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\Delta V_{tot}$ [km/s]';
ax.YLabel.FontSize = 24;
ax.YAxis.FontName= 'CMU serif';
ax.YAxis.FontSize= 16;

