function [V_0_x,V_0_y,V_0_z] = T2(Isp,g,launch_angle,m_final,m_empty)
% This function calculates the change in velocity due to the Russian
% teacher's equation. Inputs are Isp, g (m/s^2), launch_angle (degrees),
% m_final (kg), m_empty (kg)

V_0 = Isp*g*log(m_final/m_empty); % m/s
V_0_x = cosd(launch_angle) * V_0; % m/s
V_0_y = 0; % m/s [Due to the launch occuring along the x axis]
V_0_z = sind(launch_angle) * V_0; % m/s

end

