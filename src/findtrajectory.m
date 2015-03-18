% Bottle rocket code

% pressure_ambient-ambient pressure; volume_bottle-volume of bottle; discharge_coeff-discharge coefficient;
% Rhow-density of water; P_0-initial absolute pressure of air;
% 1.4-ratio of specific heats for air; g-gravity; drag_coeff-drag coefficient;
% gas_constant-gas constant for air; CL=lift coefficient of wing; V_0-initial volume
% of water; m_air_i-initial mass of air; p_end-pressure of air once water is
% expelled; WingArea-area of wing; A_b-bottle cross sectional area;
% area_throat-throat of nozzle cross sectional area; RocketMass-mass of rocket;
% T_0-initial temperature

global pressure_ambient density_h20 volume_bottle discharge_coeff P_0 gravity drag_coeff gas_constant V_0 m_air_i p_end A_b area_throat T_0
global force_x MassW0;
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

T_0 = 300.0;

m_air_i = (P_0/(gas_constant*T_0))*V_0;
MassW0 = VW0*density_h20;

% at the end of the water expulsion phase
p_end = P_0*(V_0/volume_bottle)^1.4;
%Tend = T_0*((p_end/P_0)^((1.4-1)/1.4));



EndTime = 50; % seconds
% Set up the constants


% Set up the initial conditions
vars_init(1) = 0; % initial velocity of rocket
vars_init(2) = AngI; % initial angle of attack
vars_init(3) = 0.0; % inital horizontal displacement
vars_init(4) = AltL; % initial vertical displacement
vars_init(5) = MassW0+m_air_i+MassB; % initial mass of rocket
vars_init(6) = V_0; % inital volume of air
vars_init(7) = m_air_i; % initial mass of air

vars_init=vars_init';

% Set up the integration time

time = linspace(0,5,1000);
% call ode45
[t,vars_dt] = ode45('eqns',time,vars_init); 
stopindex = 641;
stoptime = stopindex*time(2);
% temp = 9.0039 .* gravity .* log(1.1542./MassB);

figure;
plot(t,vars_dt(:,(1)));
xlabel('Time(s)');
ylabel('Velocity(m/s)');
title('Rocket Velocity as a Function of Time');

figure;
plot(t,vars_dt(:,(2)));
xlabel('Time(s)');
ylabel('Angle(rad)');
title('Rocket Angle as a Function of Time');

figure;
plot(t,vars_dt(:,(3)));
xlabel('Time(s)');
ylabel('Distance(m)');
title('Rocket Distance as a Function of Time');

figure;
plot(t,vars_dt(:,(4)));
xlabel('Time(s)');
ylabel('Height(m)');
title('Rocket Height as a Function of Time');

figure;
plot(t,vars_dt(:,(5)));
xlabel('Time(s)');
ylabel('Rocket Mass (m)');
title('Rocket Mass as a Function of Time');

figure;
plot(t,vars_dt(:,(6)));
xlabel('Time(s)');
ylabel('Air Volume');
title('Air Volume as a Function of Time');

figure;
plot(vars_dt(:,3),vars_dt(:,4),'r-')
title('Bottle Rocket Trajectory Using Validation Constants');
legend('Bottle Rocket Trajectory');
xlabel('Distance(m)');
ylabel('Height(m)');
