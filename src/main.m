%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Main Function

%% Purpose: Define constants and assess differences between models

%% Constants
global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial p_end bottle_area area_throat T_initial rho_air

%% Launch Angles and Positions
pos_initial = [0 0 0];                                                          % [m] Initial position
theta = pi/4;                                                                   % [rad] Launch angle
phi = pi/2;                                                                     % [rad] Aim heading
velocity_initial = [0 0 0];                                                     % [m/s] Initial velocity

%% Material Constants
rho_air = 1.05;                                                                 % [kg/m^3] Density of air
gravity = 9.81;                                                                 % [m/s^2] Gravitation acceleration
gas_constant = 287.15;                                                          % [Specific gas constant of air]
pressure_ambient = 101325;                                                      % [Pa] Ambient pressure
density_water = 1000.0;                                                         % [kg/m^3] Density of water
T_initial = 300.0;                                                              % [k] Ambient air temperature
velocity_initial = [0 0 0];                                                     % [m/s^2] Initial velocity
MassB = 0.15;                                                                   % [kg] Bottle mass

%% Bottle Dimensions and Weights
drag_coeff = 0.3;                                                               % [N/A] Drag coefficient of rocket
discharge_coeff = .9;                                                           % [N/A] Nozzle efficiency of the rocket
volume_bottle = 0.002;                                                          % [m^3] Volume of bottle
VW0 = 0.001;                                                                    % [m^3] Initial volume of water in bottle
volume_initial = volume_bottle-VW0;                                             % Initial volume of air in bottle

pressure_gage = 40.0;                                                           % [psi] Initial gauge pressure of air in bottle
pressure_absolute = convpres(pressure_gage,'psi','Pa');
pressure_absolute = pressure_gage+pressure_ambient;                             % [Pa] Intial absolute pressure of air in bottle

bottle_diameter=.105;                                                           % [m] Diameter of bottle
bottle_area = (pi*bottle_diameter^2);

throat_diameter=.021;                                                           % [m] Diameter of nozzle/throat
throat_area = (pi*throat_diameter^2);

mass_air_initial = (pressure_absolute/(gas_constant*T_initial))*volume_initial; % [kg] Initial mass of air
initial_mass_water = VW0*density_water;                                         % [kg] Initial mass of water

p_end = pressure_absolute*(volume_initial/volume_bottle)^1.4;                   % [Pa] Pressure at end of phase one;

mass_rocket_initial = initial_mass_water+mass_air_initial+MassB;                % [kg] initial mass of rocket



%% Numeric Integration initial constants
vars_init(1) = pos(1);
vars_init(2) = pos(2);
vars_init(3) = pos(3);
vars_init(1) = velocity_initial(1);
vars_init(1) = velocity_initial(2);
vars_init(2) = velocity_initial(3);
vars_init(5) = mass_rocket_initial;
vars_init(6) = volume_initial;
vars_init(7) = mass_air_initial;

vars_init=vars_init';

t_span = [0 10];
options = odeset('AbsTol',1e-14,'RelTol',1e-14);
[x,y,z,vx,vy,vz,mass_rocket,volume,mass_air] = ode45(@eqns,t_span,vars_init,options);

