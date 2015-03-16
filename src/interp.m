% function [m_flow, v_e, I_sp] = interp(air_density, area_throat)

% inputs: 
close all

% load static test data
loading = load('8am_3_16_test1');
load_sum = loading(:,3); % load cell summed force [lbf]
count = length(load_sum);
main_x_axis = (1:count)';

% define variables
sample_freq = 1.652*1000; % sample freq [Hz]
sample_t = 1/sample_freq; % sample interval [s]
g = 9.81; % gravitational acceleration [m/s/s]
m_fuel = 1*.001*999.97; % fuel mass [kg]

% plot static test data
figure
plot(main_x_axis,load_sum)

% identify thrust phase of sample data
index = find(load_sum>2);
index2 = find(load_sum==min(load_sum));
index3 = index(1)-2:index2;
T_data = load_sum(index3); % [lbf]
count = length(T_data);
x_axis = 1:count;

% determine accelerometer offset
ref_nom = zeros(count,1);
ref_offset = (linspace(0,min(load_sum),count))'; % [lbf]

% plot accelerometer offset with static test data
figure
plot(x_axis,T_data,'-b',x_axis,ref_nom,'-g',x_axis,ref_offset,'-r')

% adjust data for accelerometer offset and convert to N
T_data_adj = (T_data - ref_offset)*4.44822162; % [N]

% plot adjusted data
figure
plot(x_axis,T_data_adj,'-b',x_axis,ref_nom,'-g')

% determine specific impulse of adjusted data using numerical integration
x_sp = linspace(0,(count*sample_t),count);
I_sp = trapz(x_sp,T_data_adj)/(g*m_fuel); % [1/s]

% calculate propellent mass flow rate and exit velocity using thrust data
v_e = sqrt(T_data_adj/(air_density*area_throat));
m_flow = rho*A*v_e;

% end




