function [vars_dt] = eqns(t,vars, method)
%{
%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Main Function

%% Purpose: To calculate the instantaneous rate of change of each variable listed in inputs and outputs.

| Input       | Unit | Description                             |
|-------------+------+-----------------------------------------|
| position    | m    | A vector of the current rocket position |
| velocity    | m    | Current rocket velocity                 |
| volume_air  | m^3  | Air volume                              |
| mass_rocket | kg   | Rocket mass                             |
| mass_air    | kg   | Air mass                                |

| Output | Unit  | Description                   |
|--------+-------+-------------------------------|
| dposdt | m/s   | Rate of change of position    |
| dvdt   | m/s^2 | Rate of change of velocity    |
| dVdt   | m^3/s | Rate of change of volume      |
| dmrdt  | kg/s  | Rate of change of rocket mass |
| dmadt  | kg/s  | Rate of change of air mass    |


%}
    global velocity_windx velocity_windy velocity_windz pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air mass_rocket_initial launch_angle launch_rail_length mass_water_initial test_data friction_coefficient

    velocity_wind = [velocity_windx, velocity_windy, velocity_windz];
    position = [vars(1) vars(2) vars(3)];
    velocity = [vars(4) vars(5) vars(6)];
    volume_air = vars(7);
    mass_rocket   = vars(8);
    mass_air  = vars(9);

    


    switch method
      case 'thermodynamic'
        % Ignores static test stand data, water and air viscosity (friction)
        % Ignores change in bottle shape and volume due to internal pressure
        % Allows investigation of changes to basic parameters (water volume, pressure, etc.) without physical testing
        % Provides model of force_thrust force over time
        % Assumes ideal nozzle
        % Phase 1 - before all the water is expelled
        [force_thrust,dVdt,dmrdt,dmadt] = model_thermodynamic(vars);
      case 'tsiolkovsky'                                  % Second method.
                                                          % Assumes instantaneous initial acceleration/ velocity
                                                          % Higher initial force_drag due to larger velocity, less inertia during initial stages of flight
                                                          % Generally greater initial angle of trajectory due to large initial velocity (smaller initial effects of gravity)
                                                          % Simpler incorporation of statistics (use mean Isp for model)

        % No pressure, no water. Air is at ambient pressure.
        dmrdt = 0;
        dVdt = 0;
        dmadt = 0;
        force_thrust = 0;
        M_e = 0;
        v_e = 0;
        m_dot_air = 0;
        
      case 'interpolation'                                  % Third method.
                                                            % Direct use of measured time dependent force_thrust data
                                                            % More complex accounting for initial conditions (keeping rocket from falling backwards through pad)
                                                            % Uncertainty in static test stand data
                                                            % Requires static tests for any condition changes
                                                            % Assumes static performance is the same as in flight}
                                                            % No pressure, no water. Air is at ambient pressure.

        [model_thrust, m_flow] = model_interpolation(t);
        
        dmrdt = m_flow;
        dVdt = 0;
        dmadt = 0;
        force_thrust = model_thrust;
        
      otherwise
        error('None of the three cases!');
    end




    %% Trajectory equations solved for all cases

    if (position(3) <= 0)               % Rocket is on ground. Use 0 for all derivatives and return.
        for i=1:length(vars)
            vars_dt(i) = 0.0;
        end
        vars_dt = vars_dt';
        return
    else if (position(3) < launch_rail_length*sin(launch_angle)) % Rocket is still on launch rail.

            velocity_heading = [cos(launch_angle), 0, sin(launch_angle)];
            velocity_heading = velocity_heading ./ norm(velocity_heading);
            
            force_friction = mass_rocket*gravity*friction_coefficient*2;
            
            q_inf = 0.5 .* density_air * norm(velocity).^2; % Dynamic pressure
            
        else
        velocity_rel = velocity - velocity_wind;
        velocity_heading = velocity_rel ./ norm(velocity_rel); % Rocket is on neither
        force_friction = 0;
        
        q_inf = 0.5 .* density_air * norm(velocity_rel).^2; % Dynamic pressure
        
    end
    
    force_drag =  q_inf .* bottle_area .* drag_coeff; % Force_Drag due to bottle area. Only in direction of relative velocity.
    
    force_x = (force_thrust - force_drag - force_friction)*velocity_heading(1);
    force_y = (force_thrust - force_drag - force_friction)*velocity_heading(2);
    force_z = (force_thrust - force_drag - force_friction)*velocity_heading(3) + mass_rocket*gravity;
    force_vec = [force_x force_y force_z];
    dvdt = force_vec./mass_rocket;

    
    dposdt = velocity;
    

    %% Return final derivates
    vars_dt = [[dposdt]'; [dvdt]'; dVdt; dmrdt; dmadt];
    

    return
end
