function [velocity_wind] = wind_vector(magnitude1,direction1,magnitude2,direction2,height)
%% James Pavek, Noel Puldon, Haoyu Li, Jake Harrell, Nick Monahan, Aaron McCusker
%% ASEN2004 - Lab 2 - Simulated Bottle Rocket
%% Wind Vector Function

%% Purpose: Outputs a wind vector base on the inputs of magnitude of direction from sensorg

%{
| Inputs    | Units | Description                    |
|-----------+-------+--------------------------------|
| magnitude | mph   | Magnitude of wind              |
| direction | N/A   | Cardinal direction from sensor |
| 1         | N/A   | ground data                    |
| 2         | N/A   | rooftop data                   |
| height    | ft     | Instantaneous rocket height    |

| Outputs     | Units | Description                              |
|-------------+-------+------------------------------------------|
| wind_vector | m/s   | Wind vector in 3 dimensions              |
%}

%% Define constants
    theta = 20; %[deg] always CW from south as reference frame
    
    if height > 40
        magnitude = convvel(magnitude1,'mph','m/s');
        direction = direction1;
    else magnitude = convvel(magnitude2,'mph','m/s');
        direction = direction2;
    end

    %Define angle phi [deg] which starts from south and goes CW like theta
    switch direction
      case 'N'
        phi = 180;
      case 'S'
        phi = 0;
      case 'E'
        phi = 270;
      case 'W'
        phi = 90;
      case 'NW'
        phi = 225;
      case 'NE'        
        phi = 135;
      case 'SW'
        phi = 45;
      case 'SE'        
        phi = 315;
      case 'NNE'        
        phi = 202.5;
      case 'NNW'
        phi = 158.5;
      case 'SSE'
        phi = 337.5;
      case 'SSW'
        phi = 22.5;
      case 'ENE'
        phi = 248.5;
      case 'ESE'
        phi = 292.5;
      case 'WNW'
        phi = 112.5;
      case 'WSW'        
        phi = 68.5;
    end
    velocity_wind = [magnitude*-cosd(phi - theta),magnitude*sind(phi - theta), 0];

end

