function[maxDistance] = maxdistance(x);
% This file is exactly the same as findtrajectory except it returns
% the maximum distance of the params
% These params are drag coeff, launch angle, initial air pressure,
% and initial volume of water
global pressure_ambient density_h20 volume_bottle discharge_coeff P_0 gravity drag_coeff gas_constant V_0 m_air_i p_end A_b area_throat T_0
drag_coeff = x(1);
AngI = x(3);
P_0 = x(4);
VW0 = x(2);
gravity = 9.81;
gas_constant = 287.0;
discharge_coeff = .9;
pressure_ambient = 12.1*101325.0/14.7;  %8.3404e4; %100663.456;
density_h20 = 1000.0;

volume_bottle = 0.002;
V_0 = volume_bottle-VW0; %initial volume of air in bottle
VelI = 0.0;  % initial velocity of rocket
AltL = 0.01; % initial altitide of launch
MassB = 0.15; % mass of bottle itself

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



EndTime = 5; % seconds
% Set up the constants


% Set up the initial conditions
vars_init(1) = VelI; % initial velocity of rocket
vars_init(2) = AngI; % initial angle of attack
vars_init(3) = 0.0; % inital horizontal displacement
vars_init(4) = AltL; % initial vertical displacement
vars_init(5) = MassW0+m_air_i+MassB; % initial mass of rocket
vars_init(6) = V_0; % inital volume of air
vars_init(7) = m_air_i; % initial mass of air

vars_init=vars_init';

% Set up the integration time
time = linspace(0,EndTime,1000);

% call ode45
[t,vars_dt] = ode45('eqns',time,vars_init); 

maxDistance = max(vars_dt(:,3));
maxDistance = -maxDistance;
% We find the maximum distance reached with max, and then negate
% it. Because patternsolve() is a MINIMIZER, we want the lowest
% output from maxdistance, so we simply negate the output.
end
