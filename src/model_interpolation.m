function [thrust, m_flow] = model_interpolation(time)
    
%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Interpolation Model Function

%% Purpose: To use static test data to model thrust and mass flow rate over time

%{
| Inputs      | Units | Description           |
|-------------+-------+-----------------------|
| filename    | N/A   | Static test data file |
| time        | s     | Time                  |

| Outputs | Units | Description                              |
|---------+-------+------------------------------------------|
| thrust  | N     | Thrust at time in sample data            |
| m_flow  | kg/s  | Calculated mass flow rate in sample data |
%}


global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial thrust_data friction_coefficient thrust_data sample_freq

sample_t = 1/(sample_freq); % sample interval [s]

count = length(thrust_data);
x_axis = 1:count;

% determine if desired time is within thrust phase
if time > (count*sample_t)
    m_flow = 0;
    exit_velocity = 0;
    thrust = 0;
    I_sp = 0;
    return
end


x_sp = linspace(0,(count*sample_t),count);

index = (time/sample_t);
if (time ==0)
    i = 1;
    j = 2;
else 
    j = ceil(index);
    i = floor(index);
end
if (i ~= 0 && j ~=0)
    thrust = abs(thrust_data(i)+(thrust_data(j) - thrust_data(i))*((index - i)/(j-i)));
    exit_velocity = sqrt(thrust/(density_water*throat_area));
    m_flow = -density_water*throat_area*exit_velocity;
else
    thrust = 0;
    m_flow = 0;
end

        
end