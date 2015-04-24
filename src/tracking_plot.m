%Camera Data
x = [1384.38380000000;1382.48150000000;1382.48150000000;1308.29390000000;1205.57270000000;1161.82100000000;1152.30980000000;1140.89630000000;1089.53570000000;1059.09970000000;1011.54360000000;969.694200000000;956.378500000000;939.258300000000;925.942500000000;891.702100000000;884.093100000000;870.777400000000;853.657200000000;846.048200000000;821.319000000000;808.003300000000;804.198800000000;792.785300000000;779.469600000000;762.349400000000;750.935900000000;741.424700000000;728.109000000000;720.500000000000;703.379800000000;690.064100000000;682.455100000000;672.943900000000;663.432600000000;650.116900000000;638.703400000000;627.290000000000;621.583200000000;610.169700000000;596.854000000000;591.147300000000;583.538300000000;572.124800000000;560.711400000000;547.395600000000;539.786700000000;530.275400000000;520.764200000000;516.959700000000;501.741700000000;494.132800000000;490.328300000000;482.719300000000;473.208100000000;456.087800000000;444.674400000000;433.260900000000;414.238400000000;406.629500000000;400.922700000000;395.216000000000;387.607000000000;374.291300000000;370.486800000000;359.073300000000;349.562100000000;341.953100000000;326.735100000000;321.028400000000;313.419400000000;305.810400000000;298.201500000000;284.885700000000;277.276800000000;267.765500000000;258.254300000000;244.938600000000;231.622900000000;224.013900000000;214.502600000000;204.991400000000;193.577900000000;184.066700000000;176.457700000000;168.848700000000;153.630800000000;147.924000000000;132.706100000000;126.999300000000;113.683600000000;98.4657000000000;83.2477000000000;64.2252000000000;45.2028000000000;35.6915000000000;20.4736000000000;12.8646000000000;7.15790000000000];
% Horizontal_data = x./8.5 %in
datalength = length(Horizontal_data);
% rezero = ones(datalength,2);
%  rezero(:,1) = 1440.*rezero(:,1);
%  rezero(:,2) = 1080.*rezero(:,2);
y = [976.114300000000;976.114300000000;974.212000000000;919.046900000000;888.611000000000;839.152600000000;823.934600000000;812.521100000000;764.965000000000;726.920100000000;671.755000000000;626.101100000000;618.492100000000;601.371900000000;582.349400000000;550.011200000000;538.597800000000;530.988800000000;513.868600000000;502.455100000000;481.530400000000;472.019200000000;460.605700000000;452.996700000000;439.681000000000;424.463000000000;414.951800000000;403.538300000000;395.929300000000;390.222600000000;373.102400000000;367.395600000000;357.884400000000;348.373200000000;344.568700000000;333.155200000000;325.546200000000;319.839500000000;314.132800000000;308.426000000000;302.719300000000;298.914800000000;291.305800000000;293.208100000000;281.794600000000;279.892300000000;274.185600000000;274.185600000000;268.478900000000;270.381100000000;264.674400000000;266.576600000000;260.869900000000;262.772100000000;255.163100000000;255.163100000000;255.163100000000;255.163100000000;253.260900000000;257.065400000000;258.967600000000;262.772100000000;258.967600000000;260.869900000000;262.772100000000;264.674400000000;274.185600000000;276.087800000000;276.087800000000;279.892300000000;279.892300000000;285.599100000000;291.305800000000;300.817000000000;304.621500000000;314.132800000000;316.035000000000;327.448500000000;331.253000000000;336.959700000000;352.177700000000;354.079900000000;363.591100000000;373.102400000000;384.515900000000;392.124800000000;405.440600000000;413.049500000000;424.463000000000;433.974200000000;451.094500000000;472.019200000000;483.432600000000;506.259600000000;530.988800000000;542.402200000000;567.131400000000;580.447200000000;589.958400000000];
% Vertical_data = y./5.5 %in

for i = 1:datalength
   hold on  
     
    position_x = 1440 - x(i);
   % position_x1 = position_x -x(1);
    distx(i) = position_x./5.5; %ft
  
    position_y = 1080 - y(i);
    disty(i) = position_y./8.5; %ft
    
   % position_y1 = position_y - y(1);
      
end
   
 
   distx = distx - distx(1);
   disty = disty - disty(1);
   scatter(distx,disty,'*');
 
fit = polyfit(distx,disty,4);
    x_fit = linspace(0,300);
  y_fit = polyval(fit,x_fit);
    plot(x_fit,y_fit,'r');
    legend('camera data','polyfit data');
    x2 = linspace(0,300,3000);
    plot(x2,0);
    xlabel('Horizontal distance(ft)')
    ylabel('Vertical distance(ft)')
     
  hold off