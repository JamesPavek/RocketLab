function [I_sp] = model_interpolation_isp(filename)
    
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

global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data
        

        % define variables
        sample_freq = 1.652000;
        sample_t = 1/(sample_freq*1000); % sample interval [s]

        loading = load(filename);

        load_sum = loading(:,3); % load cell summed force [lbf]
                                 % count = length(load_sum);
                                 % main_x_axis = (1:count)';

        % plot static test data
        % figure
        % plot(main_x_axis,load_sum)

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
        % figure
        % plot(x_axis,T_data,'-b',x_axis,ref_nom,'-g',x_axis,ref_offset,'-r')

        % adjust data for accelerometer offset and convert to N
        T_data_adj = (T_data - ref_offset)*4.44822162; % [N]

        % plot adjusted data
        % figure
        % plot(x_axis,T_data_adj,'-b',x_axis,ref_nom,'-g')

        % determine specific impulse of adjusted data using numerical integration
        x_sp = linspace(0,(count*sample_t),count);
        I_sp = abs(trapz(x_sp,T_data_adj)./(gravity*mass_water_initial)); % [s]

        end

    