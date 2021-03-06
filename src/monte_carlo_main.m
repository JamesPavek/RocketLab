%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Monte Carlo Simulation

%% Purpose: Create a Monte Carlo Simulation of the possible landing area for the rocket

clear all
close all
clc
%% Constants

set(0,'DefaultAxesFontSize',20)

global velocity_windx velocity_windy velocity_windz pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial launch_angle launch_rail_length mass_water_initial thrust_data friction_coefficient sample_freq

tic

LoopNumber = 5000;
%% Set uncertainties
std_wind = 3; % mph
std_temp = 1; % K
std_bottle_mass = 0.005; % kg
std_launch_rail = 0; % m
std_drag = 0.075; % unitless
std_discharge = 0.1; % unitless
std_volume_water = 0.000005; % m^3
std_pressure = 2; % psi
std_launch_angle = 2 * (pi/180); % rad
std_Isp = 0.25;

%% Material and Atmospheric Constants
gravity = -9.81;                                                                          % [m/s^2] Gravitation acceleration
gas_constant = 287.15;                                                                    % [Specific gas constant of air]
pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
density_water = 1000.0;                                                                   % [kg/m^3] Density of water
temperature_initial = convtemp(60,'F','K');                                                              % [k] Ambient air temperature
velocity_initial = [0 0 0];                                                               % [m/s^2] Initial velocity
bottle_mass = 0.154;                                                                       % [kg] Bottle mass
launch_rail_length = 0.6;                                                                 % [m] Measured from cg to end of rail
friction_coefficient = 0.0;                                                               % [N/A] Friction coefficient of each launch rods (mu_k)
density_air = pressure_ambient/(gas_constant*temperature_initial);                        % [kg/m^3] Density of air

%% Bottle Dimensions and Weights
drag_coeff = 0.25;                                                                         % [N/A] Drag coefficient of rocket
discharge_coeff = .9;                                                                     % [N/A] Nozzle efficiency of the rocket
volume_bottle = 0.002;                                                                    % [m^3] Volume of bottle
volume_water_initial = 0.001;                                                             % [m^3] Initial volume of water in bottle

pressure_gage = 40.0;                                                                     % [psi] Initial gauge pressure of air in bottle
pressure_gage = convpres(pressure_gage,'psi','Pa');
pressure_absolute = pressure_gage+pressure_ambient;                                       % [Pa] Intial absolute pressure of air in bottle

bottle_diameter = .109;                                                                     % [m] Diameter of bottle
bottle_area = (pi*(bottle_diameter/2)^2);

throat_diameter=.021;                                                                     % [m] Diameter of nozzle/throat
throat_area = (pi*(throat_diameter/2)^2);

%% Thrust Data
sample_freq = 1652.00;                  % [Hz]

files = dir('data/model/Group*');
num_files = length(dir('data/model/Group*'));

for i = 1:num_files
    
    %% Read in all files in data/
    
    filename = strcat('data/model/',files(i).name);
    data_temp = load(filename, '-ascii'); % Note: for now all group tests are included
    data_temp = convforce(data_temp,'lbf','N');
    
    data(i).thrust_raw = data_temp(:,3);
    
    %% Crop to important sections (TODO)
    
    
    [~,start] = max(data(i).thrust_raw);
    start = start - 30;
    if (start < 1)
        start = 1;
    end
    
    end_data = data(i).thrust_raw(start:end);
    [~,finish] = min(end_data);
    
    finish = finish + start;
    
    count = length(start:finish);
    
    data_offset = linspace(0,min(data(i).thrust_raw),count)';
    
    data(i).thrust = data(i).thrust_raw(start:finish) - data_offset;
    
    data(i).time = (1:length(data(i).thrust))./sample_freq;
    
    data(i).total_time = (length(data(i).thrust))./sample_freq;
    
    
end

%% Mean Thrust
mean_thrust = [];

for i = 1:num_files
    for j = 1:length(data(i).thrust)
        
        value = data(i).thrust(j);
        if (i == 1)
            mean_thrust(j) = value;
        else
            
            if(length(mean_thrust)<length(data(i).thrust) && j>length(mean_thrust))
                mean_thrust(j) = 0;
            else
                mean_thrust(j) = mean_thrust(j) + value;
            end
        end
    end
    
end

thrust_data = mean_thrust ./ num_files;

thrust_data = mean_thrust ./ num_files;

final_mass = bottle_mass + mass_air_initial;

I_sp = 1.3972;

%% Thermodynamic Model Initial Conditions
pos_initial = [0 0 0.1];                                                                  % [m] Initial position
launch_angle = pi/4;                                                                      % [rad] Launch angle
velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity
mass_water_initial = volume_water_initial*density_water;                                  % [kg] Initial mass of water
volume_initial = volume_bottle-volume_water_initial;                                      % Initial volume of air in bottle
mass_air_initial = (pressure_absolute/(gas_constant*temperature_initial))*volume_initial; % [kg] Initial mass of air

pressure_end = pressure_absolute*(volume_initial/volume_bottle)^1.4;                      % [Pa] Pressure at end of phase one;
mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;                    % [kg] initial mass of rocket

t_span = [0 10];
% ode_options = odeset('AbsTol',1e-8,'RelTol',1e-8);
ode_options = odeset();

%% Interpolation Method


velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity

% Set the varibles so that the baseline isn't altered
velocity_wind1 = 0;
wind_direction = 'W';
launch_angle1 = launch_angle;
bottle_mass1 = bottle_mass;
launch_rail_length1 = launch_rail_length;
volume_water_initial1 = volume_water_initial;
temp1 = convtemp(60,'F','K');
drag_coeff1 = drag_coeff;
discharge_coeff1 = discharge_coeff;
pressure_gage1 = pressure_gage;
I_sp1 = I_sp;

for i = 1:LoopNumber
    fprintf('Loop Number: %5.0f of %5.0f \n',i,LoopNumber)
    
    % Clear Variables
    
    clear vars
    clear velocity_wind
    clear launch_angle
    clear bottle_mass
    clear launch_rail_length
    clear volume_water_initial
    clear temperature_initial
    clear drag_coeff
    clear discharge_coeff
    clear pressure_gage
    clear I_sp
    
    % Vary Wind
    
    velocity_wind2 = normrnd(velocity_wind1,std_wind);
    velocity_wind2 = wind_vector(velocity_wind2,wind_direction);
    velocity_wind2 = convvel(velocity_wind2,'mph','m/s');
    velocity_windx = velocity_wind2(1);
    velocity_windy = velocity_wind2(2);
    velocity_windz = velocity_wind2(3);
    
    % Vary Launch Angle
    
    launch_angle = normrnd(launch_angle1,std_launch_angle);
    
    % Vary Temp
    
    temperature_initial = normrnd(temp1,std_temp);
    
    % Vary Bottle Mass
    
    bottle_mass = normrnd(bottle_mass1,std_bottle_mass);
    
    % Vary Launch Rail Length
    
    launch_rail_length = normrnd(launch_rail_length1, std_launch_rail);
    
    % Vary Amount of Water
    volume_water_initial = normrnd(volume_water_initial1,std_volume_water);
    mass_water_initial = volume_water_initial * density_water;
    volume_initial = volume_bottle-volume_water_initial;
    mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;
    
    % Vary Drag
    
    drag_coeff = normrnd(drag_coeff1,std_drag);
    
    % Vary Discharge
    
    discharge_coeff = normrnd(discharge_coeff1,std_discharge);
    
    % Vary Pressure
    
    pressure_gage = normrnd(pressure_gage1,std_pressure);
    pressure_gage = convpres(pressure_gage,'psi','Pa');
    pressure_absolute = pressure_gage + pressure_ambient;
    
    % Vary Thrust Data
    
    I_sp = normrnd(I_sp1,std_Isp);
    
    vars_init(1) = pos_initial(1);
    vars_init(2) = pos_initial(2);
    vars_init(3) = pos_initial(3);
    vars_init(4) = velocity_initial(1);
    vars_init(5) = velocity_initial(2);
    vars_init(6) = velocity_initial(3);
    vars_init(7) = volume_initial;
    vars_init(8) = mass_rocket_initial;
    vars_init(9) = mass_air_initial;

    method = 'thermodynamic';
    [t,vars1] = ode45(@(t,vars1) eqns(t,vars1',method),[0 10], vars_init,ode_options);
    
    final_mass = bottle_mass + mass_air_initial;
    [vx,vy,vz] = model_tsiolkovsky(I_sp,gravity,launch_angle,mass_rocket_initial,final_mass);
    
    velocity_initial = [vx,vy,vz];                                                            % [m/s] Initial velocity
    vars_init(4) = velocity_initial(1);
    vars_init(5) = velocity_initial(2);
    vars_init(6) = velocity_initial(3);
    method = 'tsiolkovsky';
    [t,vars2] = ode45(@(t,vars2) eqns(t,vars2',method),[0 10], vars_init,ode_options);
    MaxValuesx = [vars1(end,1) vars2(end,1)];
    MaxValuesy = [vars1(end,2) vars2(end,2)];
    endpoints(i,:) = [mean(MaxValuesx) mean(MaxValuesy)];
end
endpoints;
STDx  = std(endpoints(:,1));
STDy  = std(endpoints(:,2));

Numerator = 0;
Denomx = 0;
Denomy = 0;
for i = 1:length(endpoints)
    Numerator = Numerator + ((endpoints(i,1)-mean(endpoints(:,1)))*(endpoints(i,2)-mean(endpoints(:,2))));
    Denomx = Denomx + (endpoints(i,1)-mean(endpoints(:,1)))^2;
    Denomy = Denomy + (endpoints(i,2)-mean(endpoints(:,2)))^2;
end

Correlation = Numerator/sqrt(Denomx*Denomy);

% Create Covariance Matrix
P = [STDx^2, Correlation*STDx*STDy; Correlation*STDx*STDy, STDy^2];

% Get ellipse data
h = error_ellipse(P);
x = get(h,'XData');
y = get(h,'YData');

Mean_X = mean(x);
Mean_Y = mean(y);
drift_angle = atand(mean(endpoints(:,2))/mean(endpoints(:,1)));
Average_Distance = convlength(norm(mean(endpoints(:,1)),mean(endpoints(:,2))),'m','ft');

figure
hold on
plot(x + mean(endpoints(:,1)), y + mean(endpoints(:,2)),'red');
plot(2*x + mean(endpoints(:,1)), 2*y + mean(endpoints(:,2)),'black');
plot(3*x + mean(endpoints(:,1)), 3*y + mean(endpoints(:,2)),'green');
plot(endpoints(:,1),endpoints(:,2),'b.')
title('Monte Carlo Simulation')
xlabel('Distance Down Range (ft)')
ylabel('Drift Distance (ft)');
legend('1 Standard Deviation','2 Standard Deviations','3 Standard Deviations','Landing Point');
hold off

TOC = toc;
% Print results
fprintf('\nArea of 1 std Ellipse: %4.2f ft^2',polyarea(x,y))
fprintf('\nArea of 2 std Ellipse: %4.2f ft^2',polyarea(2*x,2*y))
fprintf('\nArea of 3 std Ellipse: %4.2f ft^2',polyarea(3*x, 3*y))
fprintf('\nTotal Time to Run: %4.1f seconds \n',TOC);
fprintf('Average Time per Loop: %4.1f seconds \n',TOC/LoopNumber);