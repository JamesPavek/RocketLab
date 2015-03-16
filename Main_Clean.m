% Aaron McCusker
% Group 17
% ASEN 2004
% Lab 2
% Created and Modified February, 2015
%
% This script imports data regarding a clean F-16 in a wind tunnel. It
% takes the data and creates a drag polar based on the average data, and
% then compares that drag polar to other drag polars for the F-16.
% Data from multiple lab groups will be analyzed and combined to increase
% sample size.
% Drag polars from the software JET3 will be used, as well as data pulled
% from NASA figures regarding the drag polar for a supersonic F-16.
%
% HERRORBAR by Jos(10584) was used to create horizontal error bars. His
% code was legally obtained from 
% http://www.mathworks.com/matlabcentral/fileexchange/3963-herrorbar

close all
clear all
clc

group_count = 1;
Group_Number = 1;
scale = 1/48; % From lab document

chord = 15; % m
chord = chord*scale; % m

for GROUP = [01 02 03 04 05 06 14 15 16 17 18 19]
    Group_Num = num2str(GROUP);
    clear Wind_Data; clear No_Wind_Data; clear Net_Axial; clear Net_Pitching;
    clear q; clear Lift; clear Drag; clear Average_Lift; clear Average_Drag;
    
    % Import data as per tradition
    if GROUP > 9
        Filename_NOWIND = strcat('GROUP',Group_Num,'_CLEAN_NOWIND.csv');
        Filename_WIND = strcat('GROUP',Group_Num,'_CLEAN_WIND.csv');
    else
        Filename_NOWIND = strcat('GROUP0',Group_Num,'_CLEAN_NOWIND.csv');
        Filename_WIND = strcat('GROUP0',Group_Num,'_CLEAN_WIND.csv');
    end
    No_Wind_Data = load(Filename_NOWIND);
    Wind_Data = load(Filename_WIND);

% Find net forces (Difference between forces w/ wind and forces w/o wind)
Net_Normal = Wind_Data(:,6) - No_Wind_Data(:,6); % N
Net_Axial = Wind_Data(:,7) - No_Wind_Data(:,7); % N
Net_Pitching = Wind_Data(:,8) - No_Wind_Data(:,8); % N*m
A = 16.97; % mm
A = A/1000;
Net_Pitching = Net_Pitching - (Net_Normal*A); % N*m

% Needed values to calculate coefficients of lift and drag
Surface = 21.36; % m^2 (From Jet3, Wing Surface)
Surface = Surface*(scale^2); % m^2
q = Wind_Data(:,3); % Pa
q_Values(:,group_count) = q;

    i = 1;
    % For each group, the data must be analyzed and compiled in a usable
    % manner.
    % For each Angle of attack
    for Angle_of_Attack = -5:18
        if Angle_of_Attack < 0
            % If Angle of Attack is less than 0, the x component of the
            % normal force is parasitic of drag (ie, it subtracts from the
            % drag).
            % Also, the y component of the axial force contributes to the
            % lift (ie, it increases it).
            Lift(i,group_count) = Net_Normal(i)*cosd(abs(Angle_of_Attack)) + Net_Axial(i)*sind(abs(Angle_of_Attack));
            Drag(i,group_count) = Net_Axial(i)*cosd(abs(Angle_of_Attack)) - Net_Normal(i)*sind(abs(Angle_of_Attack));
            Moment(i,group_count) = Net_Pitching(i);
        end
        if Angle_of_Attack >= 0
            % If Angle of Attack is greater than 0, the x component of the
            % normal force contributes to drag (ie, it adds to the drag).
            % Also, the y component of the axial force subtracts from the
            % lift (ie, it reduces it).
            Lift(i,group_count) = Net_Normal(i)*cosd(Angle_of_Attack) - Net_Axial(i)*sind(Angle_of_Attack);
            Drag(i,group_count) = Net_Axial(i)*cosd(Angle_of_Attack) + Net_Normal(i)*sind(Angle_of_Attack);
            Moment(i,group_count) = Net_Pitching(i);
        end
        i = i+1;
    end

% Coefficients
Drag_Coefficient(group_count,:) = Drag(:,group_count)./(q*Surface);
Lift_Coefficient(group_count,:) = Lift(:,group_count)./(q*Surface);
Moment_Coefficient(group_count,:) = Moment(:,group_count)./(q*Surface*chord);

group_count = group_count + 1;
end

for i = 1:24
    % DETERMINING ERROR
    
    % PARTIALS
    q_term = -mean(Lift(i,:))/((mean(q_Values(i,:))^2)*Surface);
    L_term = 1/(mean(q_Values(i,:))*Surface);
    D_term = 1/(mean(q_Values(i,:))*Surface);
    M_term = 1/(mean(q_Values(i,:))*Surface*chord);
    
    % STANDARD DEVIATIONS OF MEASUREMENTS
    std_q = std(q_Values(i,:));
    std_L = std(Lift(i,:));
    std_D = std(Drag(i,:));
    std_M = std(Moment(i,:));
    
    % QUADRATURE
    Standard_Lift(i) = 2*sqrt((std_q*q_term)^2+(std_L*L_term)^2);
    Standard_Drag(i) = 2*sqrt((std_q*q_term)^2+(std_D*D_term)^2);
    Standard_Moment(i) = 2*sqrt((std_q*q_term)^2+(std_M*M_term)^2);
    
    % Finding Averages
    Clean_Drag_Coefficient_Avg(i) = mean(Drag_Coefficient(:,i));
    Clean_Lift_Coefficient_Avg(i) = mean(Lift_Coefficient(:,i));
    Clean_Moment_Coefficient_Avg(i) = mean(Moment_Coefficient(:,i));
end

% Set values of Angle of Attack
Angle_of_Attack = -5:18;

% MOMENT VS ANGLE OF ATTACK
figure
hold on
plot(Angle_of_Attack,Clean_Moment_Coefficient_Avg)
line([-5 18], [0 0],'Color','r')
hold off
title('Experimental Moment Diagram for a Clean F-16')
xlabel('Angle of Attack (Degrees)')
ylabel('C_M')

% MOMENT VS ANGLE OF ATTACK (WITH ERROR)
figure
hold on
errorbar(Angle_of_Attack,Clean_Moment_Coefficient_Avg,Standard_Moment)
line([-5 18], [0 0],'Color','r')
hold off
title('Experimental Moment Diagram for a Clean F-16')
xlabel('Angle of Attack (Degrees)')
ylabel('C_M')

% ANGLE OF ATTACK VS AVG DRAG COEFFICIENTS
figure
plot(Angle_of_Attack,Clean_Drag_Coefficient_Avg)
title('Angle of Attack vs Coefficient of Drag')
xlabel('Angle of Attack (Degrees)')
ylabel('C_D')

% ANGLE OF ATTACK VS AVG LIFT COEFFICIENTS
figure
plot(Angle_of_Attack,Clean_Lift_Coefficient_Avg)
title('Angle of Attack vs Coefficient of Lift')
xlabel('Angle of Attack (Degrees)')
ylabel('C_L')

% ANGLE OF ATTACK VS AVG LIFT COEFFICIENTS (WITH ERROR)
figure
errorbar(Angle_of_Attack,Clean_Lift_Coefficient_Avg,Standard_Lift)
xlim([-6 19])
ylim([-1 3])
xlabel('Angle of Attack (Degrees)')
ylabel('C_L')
title('Average Lift Coefficients for all Measured Angles of Attack')

% ANGLE OF ATTACK VS AVG DRAG COEFFICIENTS (WITH ERROR)
figure
errorbar(Angle_of_Attack,Clean_Drag_Coefficient_Avg,Standard_Drag)
xlim([-6 19])
ylim([0 .8])
xlabel('Angle of Attack (Degrees)')
ylabel('C_D')
title('Average Drag Coefficients for all Measured Angles of Attack')

% EXPERIMENTAL DRAG POLAR
figure
plot(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg)
title('Average Experimental Subsonic Drag Polar for a Clean F-16')
xlim([0 .8])
ylim([-.75 3])
xlabel('C_D')
ylabel('C_L')

% EXPERIMENTAL DRAG POLAR (WITH ERROR)
figure
hold on
errorbar(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg,Standard_Lift)
herrorbar(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg,Standard_Drag)
hold off
title('Average Experimental Subsonic Drag Polar for a Clean F-16')
xlim([0 .8])
ylim([-.75 3])
xlabel('C_D')
ylabel('C_L')

% ANALYTICAL SOLUTION
k1 = 0.1168;
k2 = -0.0061;
CD0 = 0.0192;
i = 1;

for Analytical_Lift_Coefficient = -.4:0.01:2
    Analytical_Drag_Coefficient(i) = CD0 + k1*(Analytical_Lift_Coefficient^2) ...
    + k2*Analytical_Lift_Coefficient;
    i = i + 1;
end
Analytical_Lift_Coefficient = -.4:0.01:2;


% ANALYITCAL DRAG POLAR COMPARISON WITH EXPERIMENTAL (WITH ERROR)
figure
hold on
H = errorbar(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg,Standard_Lift,'r');
herrorbar(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg,Standard_Drag,'r');
G = plot(Analytical_Drag_Coefficient,Analytical_Lift_Coefficient,'blue');
hold off
xlim([0 .8])
ylim([-.75 3])
title('Analytical Subsonic Drag Polar for a Clean F-16 Comparison')
xlabel('C_D')
ylabel('C_L')
legend([H G],'Experimental Solution','Analytical Solution')

% ACTUAL VALUES FROM JET3
k1_actual = 0.122;
k2_actual = -0.007;
CD0 = 0.0197;
i = 1;

for Actual_Lift_Coefficient = -.4:0.01:2
    Actual_Drag_Coefficient(i) = CD0 + k1_actual*(Actual_Lift_Coefficient^2) ...
    + k2_actual*Actual_Lift_Coefficient;
    i = i + 1;
end
Actual_Lift_Coefficient = -.4:0.01:2;

% ACTUAL VALUES DRAG POLAR COMPARED WITH EXPERIMENTAL (WITH ERROR)
figure
hold on
H = errorbar(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg,Standard_Lift,'r');
herrorbar(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg,Standard_Drag,'r');
G = plot(Actual_Drag_Coefficient,Actual_Lift_Coefficient,'blue');
hold off
xlim([0 .8])
ylim([-.75 3])
title('Actual Subsonic Drag Polar for a Clean F-16 Comparison')
xlabel('C_D')
ylabel('C_L')
legend([H G],'Experimental Solution','Actual Values')

% DATA FROM NASA FIGURES
NASA_filename = 'nasa.txt';
FID_NASA = load(NASA_filename);
NASA_Lift_Coefficient = FID_NASA(:,1);
NASA_Drag_Coefficient = FID_NASA(:,2);

% NASA DRAG POLAR
figure
plot(NASA_Drag_Coefficient,NASA_Lift_Coefficient)
% xlim([0 .6])
% ylim([-.5 2])
title('NASA Data for Supersonic Drag Polar of a Clean F-16')
xlabel('C_D')
ylabel('C_L')

% COMPARISON OF ALL DRAG POLARS
figure
hold on
plot(Clean_Drag_Coefficient_Avg,Clean_Lift_Coefficient_Avg)
plot(Actual_Drag_Coefficient,Actual_Lift_Coefficient)
plot(Analytical_Drag_Coefficient,Analytical_Lift_Coefficient)
plot(NASA_Drag_Coefficient,NASA_Lift_Coefficient)
hold off
xlim([0 .8])
ylim([-.75 3])
title('Drag Polar Comparisons for a Clean F-16')
xlabel('C_D')
ylabel('C_L')
legend('Experimental','JET3 Actual','JET3 Analytical','NASA Supersonic Data')

% CALCULATING STALL SPEED
Weight = 2.15; % N (from lab document)
V_Stall_Clean = sqrt(Weight/(max(Clean_Lift_Coefficient_Avg)*.5*1.007*Surface)) % m/s
V_Land = V_Stall_Clean*1.2
close all
Mu = 1.568*(10^(-5)); % kg/ms
Rho = 1.007; % kg/m^3
V = 25; % m/s
L = 14.32 * scale; % m (From lab doc)
Reynolds = (Rho*V*L)/Mu;
L = 14.32;
V_Actual = Mu*Reynolds/(Rho*L);
close all