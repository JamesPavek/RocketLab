% Static test report data pull and Isp calculation


global rho_air area_throat MassW0 gravity

% determine sample time
sample_freq = 1.652; % [Hz]
sample_t = 1/(sample_freq*1000); % [s]

% find .txt data files and load data
files = dir('*.txt');
for i= 1:length(files)
    data.files(i).name = load(files(i).name);
end

% pull combined load cell data, adjusts for accelerometer offsets and 
% calculates Isp for each data file using same algorithm as interp.m function
s_data = [];
for i = 1:length(files)
    data_pull = data.files(i).name(:,3);
    index = find(data_pull>2);
    index2 = find(data_pull==min(data_pull));
    index3 = index(1)-2:index2;
    T_data = data_pull(index3); % [lbf]
    count = length(T_data);
    ref_nom = zeros(count,1);
    ref_offset = (linspace(0,min(data_pull),count))';
    T_data_adj = (T_data - ref_offset)*4.44822162; % [N]
    x_sp = linspace(0,(count*sample_t),count);
    I_sp(i) = trapz(x_sp,T_data_adj)./(gravity*MassW0); % [s]
    s_data = [s_data T_data_adj];
end


