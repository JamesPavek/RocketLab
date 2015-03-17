function [vars_dt] = eqnsEDIT(time,vars)

global pressure_ambient density_h20 volume_bottle discharge_coeff P_0 gravity drag_coeff gas_constant V_0 m_air_i p_end A_b area_throat T_0 
global isp force_x

v_mag_rel = sqrt(vars(1)^2 + vars(2)^2 + vars(3)^2);
theta = vars(4);
z     = vars(8);
m_r   = vars(9);
Volume_air     = vars(10);
m_air = vars(11);
wind_x = vars(12); % m/s
wind_y = vars(13); % m/s
wind_z = vars(14); % m/s

v_x_rel = vars(1) - wind_x;
v_y_rel = vars(2) - wind_y;
v_z_rel = vars(3) - wind_z;
if v_mag_rel > 0
    heading = [v_x_rel/v_mag_rel, v_y_rel/v_mag_rel, v_z_rel/v_mag_rel];
else
    heading = [cos(theta)/(sin(theta)+cos(theta)), 0, sin(theta)/(cos(theta) + sin(theta))];
end
phi = atan(v_y_rel/v_x_rel); % [rad]

% if the rocket hits the ground stop everything
if z <= 0
    for i=1:14
        vars_dt(i) = 0.0;
        v_mag_rel = 0;
    end
    return
end

rho_air = 1.05; %density of air
dyn_p = 0.5 * rho_air * v_mag_rel^2; %dynamic pressure
D = dyn_p * A_b * drag_coeff; %drag can be changed by looking at nose cones


% Phase 1 - before all the water is expelled
if Volume_air < volume_bottle
    p = P_0*(V_0/Volume_air)^1.4;
    T = T_0*(V_0/Volume_air)^(1.4-1);
    v_e = sqrt(2*(p-pressure_ambient)/density_h20);
    m_dot_h20 = discharge_coeff*density_h20*area_throat*v_e;
    F = m_dot_h20*v_e;
    Fx = F*heading(1);
    Fy = F*heading(2);
    Fz = F*heading(3);
%     fprintf('Phase 1 \n')
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
if Volume_air >= volume_bottle
    p_air = p_end*(m_air/m_air_i)^1.4;
    rho = m_air/volume_bottle;
    T_air = p_air/(gas_constant*rho);
    p_crit = p_air*(2/(1.4+1))^(1.4/(1.4-1));
%     fprintf('Phase 2 \n')
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
    Fx = F*heading(1);
    Fy = F*heading(2);
    Fz = F*heading(3);
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
    Fx = 0;
    Fy = 0;
    Fz = 0;
    M_e = 0;
    v_e = 0;
    m_dot_air = 0;
    isp = [isp; 0];
%     fprintf('Phase 3 \n')
end

% Trajectory equations solved for all cases (4 equations)

dvxdt = (Fx - D*heading(1))/m_r;
dvydt = (Fy - D*heading(2))/m_r;
dvzdt = (Fz - D*heading(3) - m_r * gravity)/m_r;
dtheta_dt = (-gravity * cos(theta))/v_mag_rel;
dphi_dt = 0;

if v_mag_rel < 1 % stops derivatives from blowing up where rocket vel is small
    dtheta_dt = 0;
end

force_x = [force_x; F];

dxdt = v_x_rel;
dydt = v_y_rel;
dzdt = v_z_rel;

% set all forcing functions

vars_dt = [dvxdt dvydt dvzdt dtheta_dt dphi_dt dxdt dydt dzdt dmrdt dVdt dmadt 0 0 0];
vars_dt = vars_dt';
end