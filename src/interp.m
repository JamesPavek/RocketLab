clear all
close all

loading = load('8am_3_16_test1');
load_sum = loading(:,3);
count = length(load_sum);
x_axis = (1:count)';

figure
plot(x_axis,load_sum)

index = find(abs(load_sum)>0.1);
load_sum2 = load_sum(index);
count = length(load_sum2);
x_axis = 1:count;

figure
plot(x_axis,load_sum2)

value = trapz(x_axis,load_sum2);