function [thrust] = model_thermodynamic(vars_dt)

v     = vars(1);
theta = vars(2);
z     = vars(4);
m_r   = vars(5);
V     = vars(6);
m_air = vars(7);

% Phase 1 - before all the water is expelled
if V < volume_bottle
    p = P_0*(V_0/V)^1.4;
    T = T_0*(V_0/V)^(1.4-1);
    v_e = sqrt(2*(p-pressure_ambient)/density_h20);
    m_dot_h20 = discharge_coeff*density_h20*area_throat*v_e;
    thrust = m_dot_h20*v_e;
    
    % set up the 3 forcing functions
    dmrdt = -m_dot_h20;
    dVdt = discharge_coeff*area_throat*v_e;
    dmadt = 0;
    temp = ((1/gravity)*sqrt(2*abs(p-pressure_ambient)/ ...
                             density_h20));
end

% Phase 2 - after water is expelled and compressed air is expanding
if V >= volume_bottle
    p_air = p_end*(m_air/m_air_i)^1.4;
    rho = m_air/volume_bottle;
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
    m_dot_air = discharge_coeff*rho_e*area_throat*v_e;
    thrust = m_dot_air*v_e+(p_e-pressure_ambient)*area_throat;
    temp = thrust/m_dot_air;

    if (thrust > 0)
        force_x = [force_x; thrust];
    end
    if (isreal(temp))
    isp = [isp; thrust/m_dot_air];
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
    thrust = 0;
    M_e = 0;
    v_e = 0;
    m_dot_air = 0;
    isp = [isp; 0];
end
end
