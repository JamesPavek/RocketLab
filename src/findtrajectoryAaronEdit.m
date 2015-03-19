% Bottle rocket code

% pressure_ambient-ambient pressure; volume_bottle-volume of bottle; discharge_coeff-discharge coefficient;
% Rhow-density of water; P_0-initial absolute pressure of air;
% 1.4-ratio of specific heats for air; g-gravity; drag_coeff-drag coefficient;
% gas_constant-gas constant for air; CL=lift coefficient of wing; V_0-initial volume
% of water; m_air_i-initial mass of air; p_end-pressure of air once water is
% expelled; WingArea-area of wing; A_b-bottle cross sectional area;
% area_throat-throat of nozzle cross sectional area; RocketMass-mass of rocket;
% T_0-initial temperature

close all
clear all
clc

global pressure_ambient density_h20 volume_bottle discharge_coeff P_0 gravity drag_coeff gas_constant V_0 m_air_i p_end A_b area_throat T_0
global force_x;
force_x = [];
gravity = 9.81;
gas_constant = 287.0;
drag_coeff = 0.3;
discharge_coeff = .9;
pressure_ambient = 12.1*101325.0/14.7;  %8.3404e4; %100663.456;
density_h20 = 1000.0;

volume_bottle = 0.002;
VW0 = 0.001; %initial volume of water in bottle
V_0 = volume_bottle-VW0; %initial volume of air in bottle
P_0i = 40.0; %initial gauge pressure of air in bottle
P_0i = P_0i*101325.0/14.7;
P_0 = P_0i+pressure_ambient; %intial absolute pressure of air
VelI = 0.0;  % initial velocity of rocket
AngI = pi/4; % initial angle of attack
AltL = 0.01; % initial altitide of launch
MassB =0.15; % mass of bottle itself

DiaB=10.5; % diameter of bottle cm
DiaT=2.1; % diameter of nozzle/throat cm
A_b = (pi*DiaB*DiaB)*(1e-4)/4;
area_throat = (pi*DiaT*DiaT)*(1e-4)/4;

T_0 = 300.0; % K

m_air_i = (P_0/(gas_constant*T_0))*V_0;
MassW0 = VW0*density_h20;

% at the end of the water expulsion phase
p_end = P_0*(V_0/volume_bottle)^1.4;
%Tend = T_0*((p_end/P_0)^((1.4-1)/1.4));



EndTime = 50; % seconds
% Set up the constants


% Set up the initial conditions
vars_init(1) = 0.001; % initial velocity of rocket (x)
vars_init(2) = 0.001; % initial velocity of rocket (y)
vars_init(3) = 0.001; % initial velocity of rocket (z)
vars_init(4) = AngI; % Launch Angle (theta)
vars_init(5) = 0; % initial angle from launch direction (phi)
vars_init(6) = 0; % inital displacement in x direction
vars_init(7) = 0; % initial displacement in y direction
vars_init(8) = AltL; % initial displacement in z direction
vars_init(9) = MassW0+m_air_i+MassB; % initial mass of rocket
vars_init(10) = V_0; % inital volume of air
vars_init(11) = m_air_i; % initial mass of air
vars_init(12) = 2; % wind in x direction [m/s]
vars_init(13) = 3; % wind in y direction [m/s]
vars_init(14) = 0; % wind in z direction [m/s]

vars_init = vars_init';
%%
% Set up the integration time
measure_time = 5; % [sec]
frequency = 100; % [Hz]
time = linspace(0,measure_time,measure_time*frequency);
% call ode45
[t,vars_dt] = ode45('eqnsAaronEdit',time,vars_init); 
stopindex = 641;
stoptime = stopindex*time(2);
% temp = 9.0039 .* gravity .* log(1.1542./MassB);

% 
figure;
plot(t,sqrt(vars_dt(:,(1)).^2 + vars_dt(:,(2)).^2 + vars_dt(:,(3)).^2));
xlabel('Time(s)');
ylabel('Velocity(m/s)');
title('Rocket Velocity as a Function of Time');
% 
% figure;
% plot(t,vars_dt(:,(2)));
% xlabel('Time(s)');
% ylabel('Angle(rad)');
% title('Rocket Angle as a Function of Time');
% 
% figure;
% plot(t,vars_dt(:,(3)));
% xlabel('Time(s)');
% ylabel('Distance(m)');
% title('Rocket Distance as a Function of Time');
% 
% figure;
% plot(t,vars_dt(:,(8)));
% xlabel('Time(s)');
% ylabel('Height(m)');
% title('Rocket Height as a Function of Time');
% 
% figure;
% plot(t,vars_dt(:,(5)));
% xlabel('Time(s)');
% ylabel('Rocket Mass (m)');
% title('Rocket Mass as a Function of Time');
% 
% figure;
% plot(t,vars_dt(:,(6)));
% xlabel('Time(s)');
% ylabel('Air Volume');
% title('Air Volume as a Function of Time');
% 
% figure;
% plot(vars_dt(:,3),vars_dt(:,4),'r-')
% title('Bottle Rocket Trajectory Using Validation Constants');
% legend('Bottle Rocket Trajectory');
% xlabel('Distance(m)');
% ylabel('Height(m)');

% %%
clear vars_dt
clear t
m_launch = 1270; % grams
m_end = 270; % grams
Isp = 1.12;
[V_0_x, V_0_y, V_0_z] = T2(Isp,gravity,AngI,m_launch,m_end);
vars_init(1) = V_0_x;
vars_init(2) = V_0_y;
vars_init(3) = V_0_z;
vars_init(9) = m_end/1000;
vars_init(10) = volume_bottle;

[t,vars_dt] = ode45('eqnsAaronEdit',time,vars_init);
% 
% figure
% plot(t,vars_dt(:,8))
% 
% % figure
% % plot(t,vars_dt(:,4))
% 
figure
plot3(vars_dt(:,6), vars_dt(:,7), vars_dt(:,8))
xlabel('X (meters)')
ylabel('Y (meters)')
zlabel('Z (meters)')
