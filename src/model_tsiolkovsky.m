function [vx,vy,vz] = model_tsiolkovsky(isp,g,launch_angle,m_initial,m_empty)
    
%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Ideal Rocket Equation Model

%% Purpose: To calculate the initial velocity and use a ballistic trajectory.

%{
| Inputs       | Units | Description                 |
|--------------+-------+-----------------------------|
| isp          | s     | Rocket Efficiency           |
| g            | m/s^2 | Acceleration due to gravity |
| launch_angle | rad   | Launch angle                |
| m_initial    | kg    | Initial rocket mass         |
| m_empty      | kg    | Final rocket_mass           |

| Outputs | Units | Description             |
|---------+-------+-------------------------|
| vx      | m/s   | x component of velocity |
| vy      | m/s   | y component of velocity |
| vz      | m/s   | z component of velocity |
%}

    velocity_initial = abs(isp*g*log(m_initial/m_empty)); % m/s
    vx = cos(launch_angle) * velocity_initial; % m/s
    vy = 0; % m/s 
    vz = sin(launch_angle) * velocity_initial; % m/s

end

