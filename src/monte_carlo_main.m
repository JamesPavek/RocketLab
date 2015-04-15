%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Sensitivity Function

%% Purpose: Create an estimated landing area based on varying parameters

% Reset Button
clear all; close all; clc
% Begin timer
tic

% Set number of loops
Loop_Number = 4;

%% Set uncertainties
std_wind = 3; % mph
std_wind_direction = 1; % Cardinal Direction (Out of 16)
std_temp = 1; % K
std_bottle_mass = 0.001; % kg
std_launch_rail = 0.01; % m
std_drag = 0.05; % unitless
std_discharge = 0.05; % unitless
std_volume_water = 0.00001; % m^3
std_pressure = 1; % psi
std_launch_angle = 3 * (pi/180); % rad

%% Constants

global friction_coefficient thrust_data sample_freq pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data

%% Material and Atmospheric Constants
velocity_wind = wind_vector(5,'W');                                                       % [m/s] Wind vector
density_air = 1.042;                                                                      % [kg/m^3] Density of air
gravity = -9.81;                                                                          % [m/s^2] Gravitation acceleration
gas_constant = 287.15;                                                                    % [Specific gas constant of air]
pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
density_water = 1000.0;                                                                   % [kg/m^3] Density of water
temperature_initial = 300.0;                                                              % [k] Ambient air temperature
velocity_initial = [0 0 0];                                                               % [m/s^2] Initial velocity
bottle_mass = 0.15;                                                                       % [kg] Bottle mass
launch_rail_length = 0.6;                                                                 % [m] Measured from cg to end of rail
friction_coefficient = 0.0;                                                               % [N/A] Friction coefficient of each launch rods (mu_k)

%% Bottle Dimensions and Weights
drag_coeff = 0.2;                                                                         % [N/A] Drag coefficient of rocket
discharge_coeff = .9;                                                                     % [N/A] Nozzle efficiency of the rocket
volume_bottle = 0.002;                                                                    % [m^3] Volume of bottle
volume_water_initial = 0.001;                                                             % [m^3] Initial volume of water in bottle

pressure_gage = 40.0;                                                                     % [psi] Initial gauge pressure of air in bottle
pressure_gage = convpres(pressure_gage,'psi','Pa');
pressure_absolute = pressure_gage+pressure_ambient;                                       % [Pa] Intial absolute pressure of air in bottle

bottle_diameter=.109;                                                                     % [m] Diameter of bottle
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

fprintf('Thrust Data Calculated \n');

velocity_wind1 = 5;
% Direction = 'ENE';
% WindDirections = ['N  ';'NNE';'NE ';'ENE';'E  ';'ESE';'SE ';'SSE';'S  ';...
%     'SSW';'SW ';'WSW';'W  ';'WNW';'NW ';'NNW'];
% Direction_Number = find(WindDirections == Direction)
velocity_wind = 0;
for i = 1:Loop_Number
    fprintf('Loop Number: %5.0f \n',i)
    clear vars
    %% Material and Atmospheric Constants
    clear velocity_wind
    velocity_wind = velocity_wind1 + std_wind*rand(1);                                                                 % [m/hr] Wind speed
%   Direction = WindDirections(Direction_Number + std_wind_direction*randi(1))
    velocity_wind = wind_vector(velocity_wind,'E');
    velocity_wind = convvel(velocity_wind,'mph','m/s');
    density_air = 1.042;                                                                      % [kg/m^3] Density of air
    gravity = -9.81;                                                                          % [m/s^2] Gravitation acceleration
    gas_constant = 287.15;                                                                    % [Specific gas constant of air]
    pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
    density_water = 1000.0;                                                                   % [kg/m^3] Density of water
    temperature_initial = 300.0 + std_temp*rand(1);                                                              % [k] Ambient air temperature
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
    ode_options = odeset('RelTol',eps);                                                             % [m/s] Initial velocity
    t_span = [0 10];
    
    %% Interpolation Method
    
    final_mass = bottle_mass + mass_air_initial;

    I_sp = model_interpolation_isp(thrust_data);
    
    velocity_initial = [0 0 0];                                                               % [m/s] Initial velocity
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
    [t,vars] = ode45(@(t,vars) eqns(t,vars',method),t_span, vars_init,ode_options);
    
    endpoints(i,:) = [vars(end,1) vars(end,2)]; % m
    
end

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

figure
hold on
plot(endpoints(:,1),endpoints(:,2),'b*')
plot(x + Mean_X, y + Mean_Y,'red');
plot(2*x + Mean_X, 2*y + Mean_Y,'black');
plot(3*x + Mean_X, 3*y + Mean_Y,'green');
title('Monte Carlo Simulation With Error Ellipses')
xlabel('Distance Down Range (m)')
ylabel('Drift Distance');
legend('Landing Point','1 Standard of Deviation','2 Standards of Deviation','3 Standards of Deviation');
hold off

% Determine Area of Ellipses
leftpoint = [find(x == min(x)), find(y == min(y))];
rightpoint = [find(x == max(x)), find(y == min(y))];
radius1 = sqrt((rightpoint(1)-leftpoint(1))^2 + (rightpoint(2)-leftpoint(2))^2);


% % Call function to get data neccesary to plot the ellipse
% [x_ellipse,y_ellipse_1,y_ellipse_2,xrad,yrad,xmean,ymean] = drawellipse(endpoints);
%
% % Plot dots and ellipse
% figure
% hold on
% plot(endpoints(:,1),endpoints(:,2),'b.')
% plot(x_ellipse,real(y_ellipse_1),'red')
% plot(x_ellipse,real(y_ellipse_2),'red')
% xlabel('Distance Down Range (m)')
% ylabel('Distance From Initial Launch Line (m)')
% xlim([0 60]);
% ylim([0 10]);
% hold off

% Stop timer
TOC = toc; % Seconds

% Print results
fprintf('\n Area of 1 std Ellipse: %4.2f m^2 \n',pi*xrad*yrad)
fprintf('\n Area of 2 std Ellipse: %4.2f m^2 \n',pi*xrad*yrad)
fprintf('\n Area of 3 std Ellipse: %4.2f m^2 \n',pi*xrad*yrad)
fprintf('Total Time to Run: %f seconds \n',TOC);
fprintf('Average Time per Loop: %f seconds \n',TOC/Loop_Number);