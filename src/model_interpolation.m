function [thrust, m_flow] = model_interpolation(filename, sample_freq, time)

% inputs: name of data file, data sample freq. [kHz], desired time [s]
% outputs: propellant mass flow rate [kg/s], propellant exit velocity [m/s]
%          thrust [N], specific impulse of data set [s]
    
    global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind mass_rocket_initial pressure_absolute launch_angle launch_rail_length mass_water_initial test_data


    % define variables
    sample_t = 1/(sample_freq*1000); % sample interval [s]

    % load static test data
    if (isempty(test_data))
        loading = load(filename);
        test_data = loading;
    else
        loading = test_data;
    end
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

    % determine if desired time is within thrust phase
    if time > (count*sample_t)
        m_flow = 0;
        v_e = 0;
        thrust = 0;
        I_sp = 0;
        return
    end

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
    %I_sp = trapz(x_sp,T_data_adj)./(gravity*mass_water_initial); % [s]

    % calculate propellent mass flow rate and exit velocity using thrust data
    % for desired time
    index = (time/sample_t);
    if (time ==0)
        i = 1;
        j = 2;
    else 
        j = ceil(index);
        i = floor(index);
    end
    if (i ~= 0 && j ~=0)
        thrust = T_data_adj(i)+(T_data_adj(j) - T_data_adj(i))*((index - i)/(j-i));
        v_e = sqrt(thrust/(density_air*throat_area));
        m_flow = density_air*throat_area*v_e;
    else
        thrust = 0;
        m_flow = 0;
    end




    