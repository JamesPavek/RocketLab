%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Sensitivity Function

%% Purpose: Create an estimated landing area based on varying parameters

% Reset Button
close all; clear all; clc;

% Begin timer
tic

% Set number of loops
Loop_Number = 100;

%% Set uncertainties
std_wind = 1; % mph
std_temp = 1; % K
std_bottle_mass = 0.001; % kg
std_launch_rail = 0.01; % m
std_drag = 0.05; % unitless
std_discharge = 0.05; % unitless
std_volume_water = 0.000001; % m^3
std_pressure = 1; % psi
std_launch_angle = 1 * (pi/180); % rad

%% Constants

global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data

for i = 1:Loop_Number
    fprintf('Loop Number: %5.0f \n',i)
    clear vars
    %% Material and Atmospheric Constants
    velocity_wind = [5 5 0];                                                                 % [m/hr] Wind speed
    velocity_wind = convvel(velocity_wind,'mph','m/s');
    velocity_wind(1) = velocity_wind(1) + std_wind*rand(1);
    velocity_wind(2) = velocity_wind(2) + std_wind*rand(1);
    density_air = 1.042;                                                                      % [kg/m^3] Density of air
    gravity = -9.81;                                                                          % [m/s^2] Gravitation acceleration
    gas_constant = 287.15;                                                                    % [Specific gas constant of air]
    pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
    density_water = 1000.0;                                                                   % [kg/m^3] Density of water
    temperature_initial = 300.0 + std_temp*rand(1);                                                              % [k] Ambient air temperature
    velocity_initial = [0 0 0];                                                               % [m/s^2] Initial velocity
    bottle_mass = 0.15;                                                                       % [kg] Bottle mass
    bottle_mass = bottle_mass + std_bottle_mass*rand(1);
    launch_rail_length = 0.6;                                                                 % [m] Measured from cg to end of rail
    launch_rail_length = launch_rail_length + std_launch_rail*rand(1);
    
    %% Bottle Dimensions and Weights
    drag_coeff = 0.3;                                                                         % [N/A] Drag coefficient of rocket
    drag_coeff = drag_coeff + std_drag*rand(1);
    discharge_coeff = .9;                                                                     % [N/A] Nozzle efficiency of the rocket
    discharge_coeff = discharge_coeff + std_discharge*rand(1);
    volume_bottle = 0.002;                                                                    % [m^3] Volume of bottle
    volume_water_initial = 0.001;                                                             % [m^3] Initial volume of water in bottle
    volume_water_initial = volume_water_initial + std_volume_water*rand(1);
    
    pressure_gage = 40.0;                                                                     % [psi] Initial gauge pressure of air in bottle
    pressure_gage = pressure_gage + std_pressure * rand(1);
    pressure_gage = convpres(pressure_gage,'psi','Pa');
    pressure_absolute = pressure_gage+pressure_ambient;                                       % [Pa] Intial absolute pressure of air in bottle
    
    bottle_diameter = .109;                                                                     % [m] Diameter of bottle
    bottle_area = (pi*bottle_diameter^2);
    
    throat_diameter = .021;                                                                     % [m] Diameter of nozzle/throat
    throat_area = (pi*throat_diameter^2);
    
    
    
    %% Thermodynamic Model Initial Conditions
    pos_initial = [0 0 0.1];                                                                  % [m] Initial position
    launch_angle = pi/4;                                                                      % [rad] Launch angle
    launch_angle = launch_angle + std_launch_angle *rand(1);
    velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity
    mass_water_initial = volume_water_initial*density_water;                                  % [kg] Initial mass of water
    volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
    mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air
    
    pressure_end = pressure_absolute*(volume_initial/volume_bottle)^1.4;                      % [Pa] Pressure at end of phase one;
    mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;                    % [kg] initial mass of rocket
    
    t_span = [0 10]; % seconds
    ode_options = odeset('AbsTol',1e-6,'RelTol',1e-6);                                                             % [m/s] Initial velocity
    
    %% Interpolation Method
    
    
    velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity
    volume_water_initial = 0;
    mass_water_initial = volume_water_initial * density_water;                                % [kg] Initial mass of water
    volume_initial = volume_bottle;                                                           % Initial volume of air in bottle
    volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
    pressure_absolute = pressure_ambient;
    mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air
    
    pressure_end = pressure_absolute*(volume_initial/volume_bottle)^1.4;                      % [Pa] Pressure at end of phase one;
    mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;                    % [kg] initial mass of rocket
    
    vars_init(1) = pos_initial(1);
    vars_init(2) = pos_initial(2);
    vars_init(3) = pos_initial(3);
    vars_init(4) = velocity_initial(1);
    vars_init(5) = velocity_initial(2);
    vars_init(6) = velocity_initial(3);
    vars_init(7) = volume_initial;
    vars_init(8) = mass_rocket_initial;
    vars_init(9) = mass_air_initial;
    
    method = 'interpolation';
    [t,vars] = ode45(@(t,vars) eqns(t,vars',method),[0 10], vars_init,ode_options);
    
    endpoints(i,:) = [vars(end,1) vars(end,2)]; % m
    
end

% Call function to get data neccesary to plot the ellipse
[x_ellipse,y_ellipse_1,y_ellipse_2,xrad,yrad,xmean,ymean] = drawellipse(endpoints);

% Plot dots and ellipse
figure
hold on
plot(endpoints(:,1),endpoints(:,2),'b.')
plot(x_ellipse,real(y_ellipse_1),'red')
plot(x_ellipse,real(y_ellipse_2),'red')
xlabel('Distance Down Range (m)')
ylabel('Distance From Initial Launch Line (m)')
xlim([0 60]);
ylim([0 10]);
hold off

% Stop timer
TOC = toc; % Seconds

% Print results
fprintf('\n Area of Ellipse: %4.2f m^2 \n',pi*xrad*yrad)
fprintf('Total Time to Run: %f seconds \n',TOC);
fprintf('Average Time per Loop: %f seconds \n',TOC/Loop_Number);