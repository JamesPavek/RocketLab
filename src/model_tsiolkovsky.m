function [vx,vy,vz] = model_tsiolkovsky(Isp,g,launch_angle,m_final,m_empty)

    V_0 = Isp*g*log(m_final/m_empty); % m/s
    vx = cos(launch_angle) * V_0; % m/s
    vy = 0; % m/s [Due to the launch occuring along the x axis]
    vz = sin(launch_angle) * V_0; % m/s

end

