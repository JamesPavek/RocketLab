%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Sensitivity Function

%% Purpose: Demonstrate the impact of different parameters on the distance that the rocket will fly.
%           These parameters DO NOT depend on static test stand data.

%% Constants

% Reset Button
clear all; close all; clc;

% Start Stopwatch
tic

set(0,'DefaultAxesFontSize',20,'defaultlinelinewidth',1.25)

SampleNumber = 150;

maximum_distance = 400; % Longest Distance a sensitivity parameter will travel (ft)

% Open Text File to Write Data Tables to
fileID = fopen('Sensitivity_Data.txt','w');

ode_options = odeset('AbsTol',1e-8,'RelTol',1e-8);
% ode_options = odeset();

global velocity_windx velocity_windy velocity_windz pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air mass_rocket_initial launch_angle launch_rail_length mass_water_initial thrust_data friction_coefficient sample_freq

%% Material and Atmospheric Constants
velocity_wind = wind_vector(0,'W');                                                       % [m/s] Wind vector
velocity_windx = velocity_wind(1);
velocity_windy = velocity_wind(2);
velocity_windz = velocity_wind(3);
density_air = 1.042;                                                                      % [kg/m^3] Density of air
gravity = -9.81;                                                                          % [m/s^2] Gravitation acceleration
gas_constant = 287.15;                                                                    % [Specific gas constant of air]
pressure_ambient = 101325;                                                                % [Pa] Ambient pressure
density_water = 1000.0;                                                                   % [kg/m^3] Density of water
temperature_initial = 300.0;                                                              % [k] Ambient air temperature
velocity_initial = [0 0 0];                                                               % [m/s^2] Initial velocity
bottle_mass = 0.154;                                                                       % [kg] Bottle mass
launch_rail_length = 0.6;                                                                 % [m] Measured from cg to end of rail
friction_coefficient = 0.2;                                                               % [N/A] Friction coefficient of each launch rods (mu_k)

%% Bottle Dimensions and Weights
drag_coeff = 0.3;                                                                         % [N/A] Drag coefficient of rocket
discharge_coeff = .7;                                                                     % [N/A] Nozzle efficiency of the rocket
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

files = dir('data/Class Static Tests/Tested and Good/Group*');
num_files = length(dir('data/Class Static Tests/Tested and Good/Group*'));

for i = 1:num_files
    
    %% Read in all files in data/
    
    filename = strcat('data/Class Static Tests/Tested and Good/',files(i).name);
    data_temp = load(filename, '-ascii'); % Note: for now all group tests are included
    
    data_temp = convforce(data_temp,'lbf','N');
    
    data(i).thrust_raw = data_temp(:,3);
    
    %% Crop to important sections
    
    [~,start] = max(data(i).thrust_raw);
    start = start - 30;
    if (start < 1)
        start = 1;
    end
    
    end_data = data(i).thrust_raw(start:end);
    
    [~,finish] = min(end_data);
    
    if (finish < 330)                   % Short air puff causes this usally. Min/max isn't good enough to check, so check manually.
        
        marker1 = finish;
        end_data = end_data(finish+50:end);
        
        
        [~,marker2] = min(end_data);
        
        finish = marker2 + 50 + marker1;
        
    end
    
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

%% Sensitivity

% Parameters to vary: Coefficient of Drag (drag_coeff, vary from 0.15 to 0.8), Mass of rocket
% (bottle_mass, vary from 0.05 to 0.5 kg), launch angle (launch_angle, vary
% from 10 to 70 degrees), propellant density (density_water, vary from 500
% to 1500)

%% Parameter 1, Varying Bottle Mass
fprintf('\n----- Parameter 1 -----\n');
fprintf(fileID,'----- Parameter 1 -----\n');
fprintf(fileID,'Varying Bottle Mass\n');
fprintf(fileID,'\nBottle Mass (kg)     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
bottle_mass = linspace(0.1,1,SampleNumber); 
B = bottle_mass;

i = 1;

for bottle_mass = linspace(0.1,1,SampleNumber)
    clear vars1
    clear vars2
    velocity_initial = [0 0 0];
    % Bottle mass sensitivity only
    mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;
    
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
    MaxValues = [max(vars1(:,1)) max(vars2(:,1))];
    A(i) = mean(MaxValues);
    
    fprintf(fileID,'%2.2f                 %5.1f\n',B(i),A(i));
    
    clear vars
    i = i + 1;
end

A = convlength(A,'m','ft');
Data_Mass = [B', A'];

figure
hold on
plot(B(2:end),A(2:end),'b*');
ylim([0 maximum_distance])
title('Varying Bottle Mass')
xlabel('Bottle Mass (kg)');
ylabel('Distance Traveled (ft)');
hold off
clear A; clear B;

% Reset used variables to standard values to avoid messing up future calculations
bottle_mass = 0.154;
mass_rocket_initial = mass_water_initial+mass_air_initial+bottle_mass;

%% Parameter 2, Varying Launch Angle

fprintf('\n----- Parameter 2 -----\n');
fprintf(fileID,'----- Parameter 2 -----\n');
fprintf(fileID,'Varying Launch Angle\n');
fprintf(fileID,'\nLaunch Angle (Degrees)     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
launch_angle_deg = linspace(10,80,SampleNumber);
B = launch_angle_deg;

i = 1;
% Begin loop
for launch_angle_deg = linspace(10,80,SampleNumber); 
    clear vars1
    clear vars2
    % Launch angle conversion, used only for variation of launch angle
    % sensitivity analysis
    launch_angle = launch_angle_deg*(pi/180);
    
    vars_init(1) = pos_initial(1);
    vars_init(2) = pos_initial(2);
    vars_init(3) = pos_initial(3);
    vars_init(4) = velocity_initial(1);
    vars_init(5) = velocity_initial(2);
    vars_init(6) = velocity_initial(3);
    vars_init(7) = volume_initial;
    vars_init(8) = mass_rocket_initial;
    vars_init(9) = mass_air_initial;
    fprintf('Parameter 2, Loop Number: %3.0f of %3.0f \n',i,length(B));    method = 'thermodynamic';
    [t,vars1] = ode45(@(t,vars1) eqns(t,vars1',method),[0 10], vars_init,ode_options);
    
    final_mass = bottle_mass + mass_air_initial;
    [vx,vy,vz] = model_tsiolkovsky(I_sp,gravity,launch_angle,mass_rocket_initial,final_mass);
    
    velocity_initial = [vx,vy,vz];                                                            % [m/s] Initial velocity
    vars_init(4) = velocity_initial(1);
    vars_init(5) = velocity_initial(2);
    vars_init(6) = velocity_initial(3);
    method = 'tsiolkovsky';
    [t,vars2] = ode45(@(t,vars2) eqns(t,vars2',method),[0 10], vars_init,ode_options);
    MaxValues = [max(vars1(:,1)) max(vars2(:,1))];
    A(i) = mean(MaxValues);
    
    fprintf(fileID,'%2d                         %5.1f\n',B(i),A(i));
    
    clear vars
    i = i + 1;
end

A = convlength(A,'m','ft');
Data_Mass = [B', A'];

figure
hold on
plot(B(2:end),A(2:end),'b*');
ylim([0 maximum_distance])
title('Varying Launch Angle')
xlabel('Launch Angle (degrees)');
ylabel('Distance Traveled (ft)');
hold off
clear A; clear B;
% Reset used variables to standard values to avoid messing up future calculations
launch_angle = pi/4;
%% Parameter 3, Varying Drag Coefficient
fprintf('\n----- Parameter 3 -----\n');
fprintf(fileID,'----- Parameter 3 -----\n');
fprintf(fileID,'Varying Drag Coefficient\n');
fprintf(fileID,'\nDrag Coefficient     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
drag_coeff = linspace(0.1,1,SampleNumber);
B = drag_coeff;

i = 1;
% Begin loop
for drag_coeff = linspace(0.1,1,SampleNumber)
    clear vars1
    clear vars2
    
    vars_init(1) = pos_initial(1);
    vars_init(2) = pos_initial(2);
    vars_init(3) = pos_initial(3);
    vars_init(4) = velocity_initial(1);
    vars_init(5) = velocity_initial(2);
    vars_init(6) = velocity_initial(3);
    vars_init(7) = volume_initial;
    vars_init(8) = mass_rocket_initial;
    vars_init(9) = mass_air_initial;
    fprintf('Parameter 3, Loop Number: %3.0f of %3.0f \n',i,length(B));
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
    MaxValues = [max(vars1(:,1)) max(vars2(:,1))];
    A(i) = mean(MaxValues);
    fprintf(fileID,'%3.3f                %5.1f\n',B(i),A(i));
    
    clear vars
    i = i + 1;
end
A = convlength(A,'m','ft');
Data_Mass = [B', A'];

figure
hold on
plot(B(2:end),A(2:end),'b*');
ylim([0 maximum_distance])
title('Varying Drag Coefficient')
xlabel('Bottle Mass (kg)');
ylabel('Distance Traveled (ft)');
clear A; clear B;

% Reset used variables to standard values to avoid messing up future
% calculations

drag_coeff = .3;

%% Parameter 4
fprintf('\n----- Parameter 4 -----\n');
fprintf(fileID,'----- Parameter 4 -----\n');
fprintf(fileID,'Varying Propellant Density\n');
fprintf(fileID,'\nProp Density (kg/m^3)     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
density_water = linspace(200,1600,SampleNumber);
B = density_water;

i = 1;
% Begin loop. Maybe go make a sandwhich because this can take a while
for density_water = linspace(200,1600,SampleNumber)
        clear vars1
    clear vars2
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
    fprintf('Parameter 4, Loop Number: %3.0f of %3.0f \n',i,length(B));
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
    MaxValues = [max(vars1(:,1)) max(vars2(:,1))];
    A(i) = mean(MaxValues);
    fprintf(fileID,'%4.1f                          %5.1f\n',B(i),A(i));
    
    clear vars
    i = i + 1;
end

badnums = find(A > 300);
A(badnums) = 0;

A = convlength(A,'m','ft');
Data_Mass = [B', A'];

figure
hold on
plot(B(2:end),A(2:end),'b*');
ylim([0 maximum_distance])
title('Varying Propellant Density')
xlabel('Propellant Density (kg/m^3)');
ylabel('Distance Traveled (ft)');
clear A; clear B;

density_water = 1000.0;

%% Parameter 5
fprintf('\n----- Parameter 5 -----\n');
fprintf(fileID,'----- Parameter 5 -----\n');
fprintf(fileID,'Varying Propellant Volume\n');
fprintf(fileID,'\nProp Volume (m^3)     Distance Traveled (m)\n\n');
% Set parameter and fill parameter vector
volume_water_initial = linspace(0.0003,0.0018,SampleNumber);
B = volume_water_initial;

i = 1;
% Begin loop. Maybe go make a sandwhich because this can take a while
for volume_water_initial = linspace(0.0003,0.0018,SampleNumber)
    clear vars1
    clear vars2
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
    fprintf('Parameter 5, Loop Number: %3.0f of %3.0f \n',i,length(B));
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
    MaxValues = [max(vars1(:,1)) max(vars2(:,1))];
    A(i) = mean(MaxValues);
    
    fprintf(fileID,'%6.6f                          %5.1f\n',B(i),A(i));
    
    clear vars
    i = i + 1;
end

badnums = find(A > 300);
A(badnums) = 0;
A = convlength(A,'m','ft');
Data_Mass = [B', A'];

figure
hold on
plot(B(2:end),A(2:end),'b*');
ylim([0 maximum_distance])
title('Varying Propellant Volume')
xlabel('Propellant Volume (m^3)');
ylabel('Distance Traveled (ft)');
clear A; clear B;

volume_water_initial = 0.001;

%% End

fclose(fileID);

TOC = toc;
fprintf('\n Time to Completion: %3.0f minutes and %2.0f seconds \n',TOC/60,rem(TOC,60))