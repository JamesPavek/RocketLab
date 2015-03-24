function wind_vector = windvector(mag,direction)
%Description: This function outputs a wind vector in [m/s] based on the inputs of
%magnitude of direction in [mph] from sensor

%Inputs: mag = magnitude of wind ie "5" [mph], direction is one of 16
%cardinal directions N,NW,NNE, etc...
%Outputs: Wind_vector is the component form of a 1x3 array [x y z] [m/s]

%% Define constants
theta = 20; %[deg] always CW from south as reference frame
mag = mag*.44704; %[m/s] converts mph to m/s

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
wind_vector = [mag*-cosd(phi - theta) mag*sind(phi - theta) 0];

end

