function [vars_dt] = eqns(vars, method)
%{
| Input       | Unit | Description                             |
|-------------+------+-----------------------------------------|
| position    | m    | A vector of the current rocket position |
| velocity    | m    | Current rocket velocity                 |
| volume_air  | m^3  | Air volume                              |
| mass_rocket | kg   | Rocket mass                             |
| mass_air    | kg   | Air mass                                |
%}
    global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind pressure_absolute launch_angle launch_rail_length

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
        if volume_air < volume_bottle
            p = pressure_absolute*(volume_initial/volume_air)^1.4;
            T = temperature_initial*(volume_initial/volume_air)^(1.4-1);
            v_e = sqrt(2*(p-pressure_ambient)/density_water);
            m_dot_h20 = discharge_coeff*density_water*throat_area*v_e;
            force_thrust = m_dot_h20*v_e;
            
            % set up the 3 forcing functions
            dmrdt = -m_dot_h20;
            dVdt = discharge_coeff*throat_area*v_e;
            dmadt = 0;
            temp = ((1/gravity)*sqrt(2*abs(p-pressure_ambient)/ density_water));
        end

        % Phase 2 - after water is expelled and compressed air is expanding
        if volume_air >= volume_bottle
            p_air = p_end*(mass_air/mass_air_i)^1.4;
            rho = mass_air/volume_bottle;
            T_air = p_air/(gas_constant*rho);
            p_crit = p_air*(2/(1.4+1))^(1.4/(1.4-1));
            
            if p_crit > pressure_ambient
                % Flow is choked
                p_e = p_crit;
                
                M_e = 1.0;
                T_e = (2/(1.4+1))*T_air;
            else
                %Flow is not choked
                p_e = pressure_ambient;
                M_e = sqrt((2/(1.4-1))*((p_air/pressure_ambient)^((1.4-1)/1.4)-1));
                T_e = T_air*(1+(((1.4-1)/2)*(M_e)^2));
            end
            
            v_e = M_e*sqrt(1.4*gas_constant*T_e);
            rho_e = p_e/(gas_constant*T_e);
            m_dot_air = discharge_coeff*rho_e*throat_area*v_e;
            force_thrust = m_dot_air*v_e+(p_e-pressure_ambient)*throat_area;
            temp = force_thrust/m_dot_air;

            if (force_thrust > 0)
                force_x = [force_x; force_thrust];
            end
            p = p_air;
            
            % set up the 3 forcing functions
            dmrdt = -m_dot_air;
            dVdt = 0;
            dmadt = -m_dot_air;
        end

        % Phase 3 - Ballistic motion
        if p <= (pressure_ambient)
            dmrdt = 0;
            dVdt = 0;
            dmadt = 0;
            force_thrust = 0;
            M_e = 0;
            v_e = 0;
            m_dot_air = 0;
        end
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
            
                q_inf = 0.5 .* density_air * norm(velocity).^2; % Dynamic pressure
            
        else
            velocity_rel = velocity - velocity_wind;
            velocity_heading = velocity_rel ./ norm(velocity_rel); % Rocket is on neither
        
            q_inf = 0.5 .* density_air * norm(velocity_rel).^2; % Dynamic pressure
        
        end
        force_drag =  q_inf .* bottle_area .* drag_coeff; % Force_Drag due to bottle area. Only in direction of relative velocity.
        force_x = (force_thrust-force_drag)*velocity_heading(1);
        force_y = (force_thrust-force_drag)*velocity_heading(2);
        force_z = (force_thrust-force_drag)*velocity_heading(3) + mass_rocket*gravity;
        force_vec = [force_x force_y force_z];
        dvdt = force_vec./mass_rocket;

    
        dposdt = velocity;
   

        %% Return final derivates
        vars_dt = [[dposdt]'; [dvdt]'; dVdt; dmrdt; dmadt];

        return
    end
