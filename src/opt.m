pressure_ambient = 86184;     
burst_pressure = 400000;
lower_bound = [0.3, 0.0001,  0, pressure_ambient]
upper_bound = [0.4, 0.002,  pi/2, burst_pressure]

% Patternsearch is a minimizer that finds the lowest value between
% these bounds. The first parameter is the function to call, the
% second is the initial guess(which are general parameters) and the
% last two are the upper and lower bounds of the vars.
max_height_vars = patternsearch(@maxheight,[0.3,0.001,pi/4,367000],[],[],[],[],lower_bound,upper_bound)
max_distance_vars = patternsearch(@maxdistance,[0.3,0.001,pi/4,367000],[],[],[],[],lower_bound,upper_bound)
disp('Vars for Max Height: Drag Coeff, Water Volume, Launch Angle, Bottle Pressure');
disp(max_height_vars);
disp('Vars for Max Distance: Drag Coeff, Water Volume, Launch Angle, Bottle Pressure');
disp(max_distance_vars);
