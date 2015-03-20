function [x,y1,y2,a,b,h,k] = drawellipse(values)
% This function generates the data neccessary to draw an ellipse around a
% set of given values.

if length(values) <= 1
    error('Not enough values to create an ellipse');
end

x_values = values(:,1);
y_values = values(:,2);

h = mean(x_values);
k = mean(y_values);

a = (max(x_values) - min(x_values))/2;
b = (max(y_values) - min(y_values))/2;

x = linspace(min(x_values),max(x_values),100);

y1 = b*sqrt(1 - ((x-h).^2)/(a^2)) + k;
y2 = -b*sqrt(1 - ((x-h).^2)/(a^2)) + k;

end

