function [vars_dt] = eqns(t_span,vars)

global pressure_ambient density_h20 volume_bottle discharge_coeff P_0 gravity drag_coeff gas_constant V_0 m_air_i p_end A_b area_throat T_0 
global isp force_x rho_air
v     = vars(1);
theta = vars(2);
z     = vars(4);
m_r   = vars(5);
V     = vars(6);
m_air = vars(7);

% if the rocket hits the ground stop everything
if z <= 0
    for i=1:7
        vars_dt(i) = 0.0;
        v=0;
    end
    return
end

rho_air = 1.05; %density of air
dyn_p = 0.5 * rho_air * v.^2; %dynamic pressure
D = dyn_p * A_b * drag_coeff; %drag can be changed by looking at nose cones


% Phase 1 - before all the water is expelled
if V < volume_bottle
    p = P_0*(V_0/V)^1.4;
    T = T_0*(V_0/V)^(1.4-1);
    v_e = sqrt(2*(p-pressure_ambient)/density_h20);
    m_dot_h20 = discharge_coeff*density_h20*area_throat*v_e;
    F = m_dot_h20*v_e;
    
    % set up the 3 forcing functions
    dmrdt = -m_dot_h20;
    dVdt = discharge_coeff*area_throat*v_e;
    dmadt = 0;
    temp = ((1/gravity)*sqrt(2*abs(p-pressure_ambient)/ ...
                             density_h20));
    if (isreal(temp))
    isp = [isp; ((1/gravity)*sqrt(2*abs(p-pressure_ambient)/density_h20))];
    end
    if (F > 0)
        force_x = [force_x; F];
    end
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
    F = m_dot_air*v_e+(p_e-pressure_ambient)*area_throat;
    temp = F/m_dot_air;
    if (F > 0)
        force_x = [force_x; F];
    end
    if (isreal(temp))
    isp = [isp; F/m_dot_air];
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
    F = 0;
    M_e = 0;
    v_e = 0;
    m_dot_air = 0;
    isp = [isp; 0];
end

% Trajectory equations solved for all cases (4 equations)
dvdt = (F - D - m_r * gravity * sin(theta))/m_r;
dtheta_dt = (-gravity * cos(theta))/v;
if v < 1 % stops derivatives from blowing up where rocket vel is small
    dtheta_dt = 0;
end
force_x = [force_x; F];
dxdt = v * cos(theta);
dzdt = v * sin(theta);
% set all forcing functions
vars_dt = [dvdt dtheta_dt dxdt dzdt dmrdt dVdt dmadt];
vars_dt = vars_dt';
end