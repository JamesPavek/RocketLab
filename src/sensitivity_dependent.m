%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Sensitivity Function

%% Purpose: Determine Sensitivity of Paramters that do depend on static test stand data

% Reset Button
close all; clear all; clc;

% Start Stopwatch
tic

% Open Text File to Write Data Tables to
fileID = fopen('Dependent_Sensitivity_Data.txt','w');

%% Constants

global Test_Stand_File pressure_ambient density_water volume_bottle pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data

t_span = [0 10];
ode_options = odeset('AbsTol',1e-7,'RelTol',1e-7);

%% Material and Atmospheric Constants
wind_speed = 0;
wind_direction = 'NW';
velocity_wind = windvector(wind_speed,wind_direction);
velocity_wind = convvel(velocity_wind,'mph','m/s');
density_air = 1.042;                                                                      % [kg/m^3] Density of air
gravity = -9.81;                                                                          % [m/s^2] Gravitation acceleration
gas_constant = 287.15;                                                                    % [Specific gas constant of air]
pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
density_water = 1000.0;                                                                   % [kg/m^3] Density of water
temperature_initial = 300.0;                                                              % [k] Ambient air temperature
velocity_initial = [0 0 0];                                                               % [m/s^2] Initial velocity
bottle_mass = 0.15;                                                                       % [kg] Bottle mass
launch_rail_length = 0.6;                                                                 % [m] Measured from cg to end of rail
pos_initial = [0 0 0.1];                                                                  % [m] Initial position
launch_angle = pi/4;                                                                      % [rad] Launch angle

%% Bottle Dimensions and Weights
drag_coeff = 0.3;                                                                         % [N/A] Drag coefficient of rocket
discharge_coeff = .9;                                                                     % [N/A] Nozzle efficiency of the rocket
volume_bottle = 0.002;                                                                    % [m^3] Volume of bottle
volume_water_initial = 0.00101;                                                             % [m^3] Initial volume of water in bottle

pressure_gage = 40.0;                                                                     % [psi] Initial gauge pressure of air in bottle
pressure_gage = convpres(pressure_gage,'psi','Pa');
pressure_absolute = pressure_gage+pressure_ambient;                                       % [Pa] Intial absolute pressure of air in bottle

bottle_diameter=.109;                                                                     % [m] Diameter of bottle
bottle_area = (pi*bottle_diameter^2);

throat_diameter=.021;                                                                     % [m] Diameter of nozzle/throat
throat_area = (pi*throat_diameter^2);

volume_water_initial = 0;
mass_water_initial = volume_water_initial * density_water;                                % [kg] Initial mass of water
volume_initial = volume_bottle;                                                           % Initial volume of air in bottle
volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
pressure_absolute = pressure_ambient;
mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air

pressure_end = pressure_absolute*(volume_initial/volume_bottle)^1.4;                      % [Pa] Pressure at end of phase one;
mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;                    % [kg] initial mass of rocket

%% Sensitivity

% Parameters to vary: Coefficient of Drag (drag_coeff, vary from 0.15 to 0.8), Mass of rocket
% (bottle_mass, vary from 0.05 to 0.5 kg), launch angle (launch_angle, vary
% from 10 to 70 degrees), propellant density (density_water, vary from 500
% to 1500)

%% Parameter 1, Varying Water Volume
fprintf('\n----- Parameter 1 -----\n');
fprintf(fileID,'----- Parameter 1 -----\n');
fprintf(fileID,'Varying Water Volume\n');
fprintf(fileID,'\nVolume Water (ml)     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
volume_water_initial = [750 1000];
B = volume_water_initial;

figure

i = 1;

for volume_water_initial = [750 1000];

% Assign file name to be read (Please use my naming convention)
Test_Stand_File = strcat('Static_Volume_Water_',num2str(volume_water_initial));
volume_water_initial = volume_water_initial / (10^6);

mass_water_initial = volume_water_initial * density_water;                                % [kg] Initial mass of water
volume_initial = volume_bottle;                                                           % Initial volume of air in bottle
volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;
mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air

vars_init(1) = pos_initial(1);
vars_init(2) = pos_initial(2);
vars_init(3) = pos_initial(3);
vars_init(4) = velocity_initial(1);
vars_init(5) = velocity_initial(2);
vars_init(6) = velocity_initial(3);
vars_init(7) = volume_initial;
vars_init(8) = mass_rocket_initial;
vars_init(9) = mass_air_initial;
fprintf('Parameter 1, Loop Number: %3.0f of %3.0f \n',i,length(B));
method = 'interpolation_sensitivity';
[t,vars] = ode45(@(t,vars) eqns(t,vars',method),[0 10], vars_init,ode_options);

A(i) = max(vars(:,1));

fprintf(fileID,'%5d                 %5.1f\n',B(i),A(i));

clear vars
i = i + 1;
end

plot(B,A,'b*');
title('Water Volume')
xlabel('Water Volume (ml)');
ylabel('Distance Traveled (m)');
xlim([(min(B) - 200) (max(B) + 200)])
clear A; clear B;

% Reset used variables to standard values to avoid messing up future
% calculations
volume_water_initial = 0.001;

%% Parameter 2, Varying Propellent Density
fprintf('\n----- Parameter 2 -----\n');
fprintf(fileID,'----- Parameter 2 -----\n');
fprintf(fileID,'Varying Propellant Density\n');
fprintf(fileID,'\nDensity Water (kg/m^3)     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
density_water = [1000];
B = density_water;

figure

i = 1;

for density_water = [1000];

% Assign file name to be read (Please use my naming convention)
Test_Stand_File = strcat('Static_Density_Water_',num2str(density_water));

mass_water_initial = volume_water_initial * density_water;                                % [kg] Initial mass of water
volume_initial = volume_bottle;                                                           % Initial volume of air in bottle
volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;
mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air

vars_init(1) = pos_initial(1);
vars_init(2) = pos_initial(2);
vars_init(3) = pos_initial(3);
vars_init(4) = velocity_initial(1);
vars_init(5) = velocity_initial(2);
vars_init(6) = velocity_initial(3);
vars_init(7) = volume_initial;
vars_init(8) = mass_rocket_initial;
vars_init(9) = mass_air_initial;
fprintf('Parameter 1, Loop Number: %3.0f of %3.0f \n',i,length(B));
method = 'interpolation_sensitivity';
[t,vars] = ode45(@(t,vars) eqns(t,vars',method),[0 10], vars_init,ode_options);

A(i) = max(vars(:,1));

fprintf(fileID,'%5d                      %5.1f\n',B(i),A(i));

clear vars
i = i + 1;
end

plot(B,A,'b*');
title('Varying Propellant Density')
xlabel('Propellant Density (kg/m^3)');
ylabel('Distance Traveled (m)');
clear A; clear B;

% Reset used variables to standard values to avoid messing up future
% calculations
density_water = 1000;

%% End
TOC = toc;
fprintf('\n Time to Completion: %3.0f minutes and %2.0f seconds \n',TOC/60,rem(TOC,60))