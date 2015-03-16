function [vars_dt] = eqns(vars, method)

global pressure_ambient density_water volume_bottle discharge_coeff P_0 gravity drag_coeff gas_constant V_0 m_air_i p_end A_b area_throat T_0 rho_air
direction = [vars(1) vars(2) vars(3)];
v = [vars(4) vars(5) vars(6)];
theta = vars(2);
phi = vars(3);
v = vars[vars(6),vars(7),vars(8)];
m_r   = vars(5);

% if the rocket hits the ground stop everything
if z <= 0
    for i=1:length(i)
        vars_dt(i) = 0.0;
        v=0;
    end
    return
end

if (method == 1)
    thrust = model_thermodynamic(vars);
else if (method == 2)
    delta_v = method2
else if (method == 3)
      
end

dyn_p = 0.5 * rho_air * v.^2; %dynamic pressure
D = dyn_p * A_b * drag_coeff; %drag can be changed by looking at nose cones


% Trajectory equations solved for all cases (4 equations)
dvdt = (F - D - m_r .* gravity)/m_r;
dtheta_dt = (-gravity * cos(theta))/v;

if v < 1 % stops derivatives from blowing up where rocket vel is small
    dtheta_dt = 0;
end


dxdt = v * cos(theta);
dzdt = v * sin(theta);

% set all forcing functions
vars_dt = [dvdt dtheta_dt dxdt dzdt dmrdt dVdt dmadt];
vars_dt = vars_dt';
end
