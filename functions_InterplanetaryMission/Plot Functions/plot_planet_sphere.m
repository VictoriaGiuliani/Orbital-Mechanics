function plot_planet_sphere(mjd2000,id,R)
% PLOT PLANET SPHERE: plots the selected planet as a 3D sphere, in the case
%   of Mercury, Venus, Earth and the Sun a texture will be added to the
%   surface.
%   NOTE: if mjd2000 is empty the function will assume that the user means
%   to plot the planetocentric view of the planet: unit of measure will
%   pass from AU to km, and the planet will be plotted at [0,0,0].
%__________________________________________________________________________   
% PROTOTYPE:
%    plot_planet_sphere(time,id,R)
% or plot_planet_sphere(time,id)
% 
% INPUT:
%   mjd2000[1]  time in MJD2000, specifies where in the solar system
%               the sphere will have its centre.                    [days]
%               If it's empty, then the center of the sphere will be
%               at [0,0,0]
%   id[1]       planet identifier, 0 for Sun                        [-]
%   R[1]        planet radius, useful to make it bigger in the solar
%               system plot                                         [km]
%
% OUTPUT:
%   figure containing the 3D sphere
%__________________________________________________________________________ 
% CONTRIBUTORS:
%   Victoria Katia Giuliani     Deepika Sampath Kumar          
%   Alberto Giuseppe Lunghi     Giulio Pelenghi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2 && id ~= 0
    R = astroConstants(id+20);
    
elseif nargin ~= 3
   error('Wrong number of input arguments') 
end


if isempty(mjd2000) || id == 0
    r = [0,0,0];
else
    [kep, ksun] = uplanet(mjd2000, id);
    [r,~] = kep2car(kep, ksun);
end

switch id
    case 0
        if nargin == 2
            R = astroConstants(3);
        end
        texture = '2k_sun.jpg';
    case 1
        texture = '2k_mercury.jpg';
    case 2
        texture = '2k_venus_atmosphere.jpg';        
    case 3
        texture = '2k_earth_daymap.jpg';
end

if ~isempty(mjd2000)
    R = R/astroConstants(2);
end
r = r/astroConstants(2);

[X,Y,Z] = sphere(50);
X= X*R + r(1); 
Y=-Y*R + r(2); % account for the flip in Y during warp
Z= Z*R + r(3);
map=imread(texture);
warp(X,Y,Z,map);
set(gca,'Ydir','normal','Xdir','normal');
    
    
    
    
    
    
    
    
    
    
    
    
    