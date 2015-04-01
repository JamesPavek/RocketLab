function [height] = apex(S1, S2, S3, alpha, d)

% inputs: angle [deg] and shoulder height [in] for positions S1, S2, S3,
%         drift angle [deg] and flight distance [ft]
%
% outputs: trajectory mean apex height [m]

ft2m = 0.3048; % ft to meters conversion factor

% convert shoulder heights to [m]
S1(2) = S1(2)/12*ft2m;
S2(2) = S2(2)/12*ft2m;
S3(2) = S3(2)/12*ft2m;

% determine distance to apex from launch pad [m]
y_vector = [0 d/2]*ft2m;

% establish x vectors for each position while accounting for drift angle [m]
x_vector = [-40 d/2*sind(alpha); 40 d/2*sind(alpha); 80 d/2*sind(alpha)]*ft2m;

% determine distance from position to apex in [m]
S1_dist = hypot(x_vector(1,1)-x_vector(1,2),y_vector(2)*cosd(alpha));
S2_dist = hypot(x_vector(2,1)-x_vector(2,2),y_vector(2)*cosd(alpha));
S3_dist = hypot(x_vector(3,1)-x_vector(3,2),y_vector(2)*cosd(alpha));

% determine apex height for each position [m]
S1_height = S1_dist*tand(S1(1))+S1(2);
S2_height = S2_dist*tand(S2(1))+S2(2);
S3_height = S3_dist*tand(S3(1))+S3(2);

% establish z vector for each position [m]
S1_zvector = [S1(2) S1_height];
S2_zvector = [S2(2) S2_height];
S3_zvector = [S3(2) S3_height];

% calculate and return mean apex height from all three positions [m]
height = mean([S1_height S2_height S3_height]);

% determine apex position [m]
x_pos = d*ft2m*sind(alpha)/2;
y_pos = d*ft2m*cosd(alpha)/2;

% plot apex position vectors
h = plot3(x_vector(1,:),y_vector,S1_zvector,'-b',x_vector(2,:),y_vector,S2_zvector,'-g',x_vector(3,:),y_vector,S3_zvector,'-c',x_pos,y_pos,height,'ro');
grid on
set(h(4),'MarkerEdgeColor','none','MarkerFaceColor','r')
xlabel('Baseline Distance [m]');
ylabel('Downrange Distance [m]');
zlabel('Height [m]');
title('Data Sheet Apex Determination Plot');
legend('S1','S2','S3','Apex','Location','SouthOutside');



end