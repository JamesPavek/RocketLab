%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Main Function

%% Purpose: Define constants and assess differences between models

%% Constants

global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data


%% Material and Atmospheric Constants
velocity_wind = [1 1 0];                                                                  % [m/hr] Wind speed
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

%% Bottle Dimensions and Weights
drag_coeff = 0.3;                                                                         % [N/A] Drag coefficient of rocket
discharge_coeff = .9;                                                                     % [N/A] Nozzle efficiency of the rocket
volume_bottle = 0.002;                                                                    % [m^3] Volume of bottle
volume_water_initial = 0.001;                                                             % [m^3] Initial volume of water in bottle

pressure_gage = 40.0;                                                                     % [psi] Initial gauge pressure of air in bottle
pressure_gage = convpres(pressure_gage,'psi','Pa');
pressure_absolute = pressure_gage+pressure_ambient;                                       % [Pa] Intial absolute pressure of air in bottle

bottle_diameter=.109;                                                                     % [m] Diameter of bottle
bottle_area = (pi*bottle_diameter^2);

throat_diameter=.021;                                                                     % [m] Diameter of nozzle/throat
throat_area = (pi*throat_diameter^2);



%% Thermodynamic Model Initial Conditions
pos_initial = [0 0 0.1];                                                                  % [m] Initial position
launch_angle = pi/4;                                                                      % [rad] Launch angle
velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity
mass_water_initial = volume_water_initial*density_water;                                  % [kg] Initial mass of water
volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
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

t_span = [0 10];
ode_options = odeset('AbsTol',1e-14,'RelTol',1e-14);

method = 'thermodynamic';
[t,vars] = ode45(@(t,vars) eqns(t,vars',method),[0 10], vars_init,ode_options);

figure; 
hold on;
set(gca,'DefaultTextInterpreter', 'latex');
set(gca,'fontsize',18);
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');

title('Rocket Trajectory'); 
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

plot3(vars(:,1),vars(:,2),vars(:,3));





%% Tsiolkovsky Case

% Find initial velocity based on the calculated ISP and ideal rocket equation.
final_mass = bottle_mass;

[vx,vy,vz] = model_tsiolkovsky(1.12,gravity,launch_angle,mass_rocket_initial,final_mass);

velocity_initial = [vx vy vz];                                                            % [m/s] Initial velocity
volume_water_initial = 0;
mass_water_initial = 0;
volume_initial = volume_bottle;                                                           % Initial volume of air in bottle
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

method = 'tsiolkovsky';
[t,vars] = ode45(@(t,vars) eqns(t,vars',method),[0 10], vars_init,ode_options);

plot3(vars(:,1),vars(:,2),vars(:,3));

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

plot3(vars(:,1),vars(:,2),vars(:,3));

legend('Thermodynamic','Tsiolkovsky','Interpolation');