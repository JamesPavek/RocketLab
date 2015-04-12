function [I_sp] = model_interpolation_isp(data)
    
%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Interpolation Model Function

%% Purpose: To use static test data to calculate Isp

%{
| Inputs      | Units | Description           |
|-------------+-------+-----------------------|
| filename    | N/A   | Static test data file |

| Outputs | Units | Description              |
|---------+-------+--------------------------|
| isp     | s     | Rocket thrust efficiency |
%}

    global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data sample_freq
        

    sample_t = 1/(sample_freq); % sample interval [s]


    % identify thrust phase of sample data
    count = length(data);
    x_axis = 1:count;


    % determine specific impulse of adjusted data using numerical integration
    x_sp = linspace(0,(count*sample_t),count);
    I_sp = abs(trapz(x_sp,data)./(gravity*mass_water_initial)); % [s]

end

