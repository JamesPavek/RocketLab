function [force_thrust,dVdt,dmrdt,dmadt] = model_thermodynamic(vars)
% Ignores static test stand data, water and air viscosity (friction)
% Ignores change in bottle shape and volume due to internal pressure
% Allows investigation of changes to basic parameters (water volume, pressure, etc.) without physical testing
% Provides model of force_thrust force over time
% Assumes ideal nozzle
% Phase 1 - before all the water is expelled
    global pressure_ambient density_water volume_bottle discharge_coeff pressure_absolute gravity drag_coeff gas_constant volume_initial mass_air_initial pressure_end bottle_area throat_area temperature_initial density_air velocity_wind pressure_absolute launch_angle launch_rail_length
    position = [vars(1) vars(2) vars(3)];
    velocity = [vars(4) vars(5) vars(6)];
    volume_air = vars(7);
    mass_rocket   = vars(8);
    mass_air  = vars(9);

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
        p_air = pressure_end*(mass_air/mass_air_initial)^1.4;
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
    
end
