%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Main Function

%% Purpose: Define constants and assess differences between models

%% Constants

global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length


%% Material and Atmospheric Constants
velocity_wind = [1 1 0];                                                                  % [m/hr] Wind speed
density_air = 1.05;                                                                       % [kg/m^3] Density of air
gravity = -9.81;                                                                           % [m/s^2] Gravitation acceleration
gas_constant = 287.15;                                                                    % [Specific gas constant of air]
pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
density_water = 1000.0;                                                                   % [kg/m^3] Density of water
temperature_initial = 300.0;                                                              % [k] Ambient air temperature
velocity_initial = [0 0 0];                                                               % [m/s^2] Initial velocity
bottle_mass = 0.15;                                                                       % [kg] Bottle mass
launch_rail_length = 0.6;               % [m] Measured from cg to end of rail

%% Bottle Dimensions and Weights
drag_coeff = 0.3;                                                                         % [N/A] Drag coefficient of rocket
discharge_coeff = .9;                                                                     % [N/A] Nozzle efficiency of the rocket
volume_bottle = 0.002;                                                                    % [m^3] Volume of bottle
VW0 = 0.001;                                                                              % [m^3] Initial volume of water in bottle

pressure_gage = 40.0;                                                                     % [psi] Initial gauge pressure of air in bottle
pressure_absolute = convpres(pressure_gage,'psi','Pa');
pressure_absolute = pressure_gage+pressure_ambient;                                       % [Pa] Intial absolute pressure of air in bottle

bottle_diameter=.109;                                                                     % [m] Diameter of bottle
bottle_area = (pi*bottle_diameter^2);

throat_diameter=.021;                                                                     % [m] Diameter of nozzle/throat
throat_area = (pi*throat_diameter^2);


%% Initial Conditions
pos_initial = [0 0 0.1];                                                                    % [m] Initial position
launch_angle = pi/4;                                                                             % [rad] Launch angle
velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity
mass_water_initial = VW0*density_water;                                                   % [kg] Initial mass of water
volume_initial = volume_bottle-VW0;                                                       % Initial volume of air in bottle
mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air

pressure_end = pressure_absolute*(volume_initial/volume_bottle)^1.4;                      % [Pa] Pressure at end of phase one;
mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;                    % [kg] initial mass of rocket

%% Numeric Integration initial constants
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

%method = 'thermodynamic';
%[t,vars] = ode45(@(t,vars) eqns(vars',method),[0 10], vars_init,ode_options);

%% Tsiolkovsky Case

[vx,vy,vz] = model_tsiolkovsky(1.12,gravity,launch_angle,mass_rocket_initial,bottle_mass);

pos_initial = [0 0 0.1];                                                                    % [m] Initial position
launch_angle = pi/4;                                                                             % [rad] Launch angle

velocity_initial = [vx vy vz];                                                               % [m/s] Initial velocity
mass_water_initial = VW0 * density_water;                                                   % [kg] Initial mass of water
volume_initial = volume_bottle;                                                       % Initial volume of air in bottle
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
[t,vars] = ode45(@(t,vars) eqns(vars',method),[0 10], vars_init,ode_options);
xpos = vars(:,1);
ypos = vars(:,2);
zpos = vars(:,3);
plot3(xpos,ypos,zpos);
