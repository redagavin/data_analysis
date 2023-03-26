clc;
clear all;

%% loading and parameter
load Spike_Data_Separated.mat
load Position_Data_Separated.mat

bin_size = 1;
stride = 0.1;

%% Bin
Spike_Data_Bin = {};
Postion_Data_Bin = {};
Speed_Data_Bin = {};

% for each run-through, we bin the data
for k = 1:length(Position_Data_Sep)

    % get the data
    Spike_Data = Spike_Data_Sep{k};
    Position_Data = Position_Data_Sep{k};

    % define the center time of the first bin and 
    % the center time of the last bin
    start_time = Spike_Data(1, 1) + bin_size/2;
    end_time = Spike_Data(end, 1) - bin_size/2;
    
    % construct the center times of all bins
    bin_center = start_time:stride:end_time;
    bin_center = bin_center';
    
    % pre-fill the arrays with zeros
    % they will be filled with real data later in the loop
    Spikes = zeros(length(bin_center), max(unique(Spike_Data(:, 2))));
    Position = zeros(length(bin_center), 1);
    Speed = zeros(length(bin_center), 1);
    
    for i = 1:length(bin_center)
    
        % get the upper bound and lower bound of the bin
        upper_bound = bin_center(i) - bin_size/2;
        lower_bound = bin_center(i) + bin_size/2;
    
        % get all the data inside the bin
        bin_spike = Spike_Data( ...
            Spike_Data(:,1)>=upper_bound & Spike_Data(:,1)<=lower_bound, :);
        bin_position = Position_Data( ...
            Position_Data(:,1)>=upper_bound & Position_Data(:,1)<=lower_bound, [1, 2]);
        bin_speed = Position_Data( ...
            Position_Data(:,1)>=upper_bound & Position_Data(:,1)<=lower_bound, [1, 5]);
        
        % process the spike data in the bin
        for j = bin_spike(:, 2)'
            Spikes(i, j) = Spikes(i, j) + 1;
        end
    
        % process the position data in the bin
        last_pi = find(bin_position(:, 1)<bin_center(i), 1, "last");
        next_pi = find(bin_position(:, 1)>bin_center(i), 1);
        position_nearby = bin_position([last_pi, next_pi], :);
        pos_lm = fitlm(position_nearby(:, 1), position_nearby(:, 2));
        center_pos = predict(pos_lm, bin_center(i));
        Position(i) = center_pos;
    
        % process the speed data in the bin
        last_si = find(bin_speed(:, 1)<bin_center(i), 1, "last");
        next_si = find(bin_speed(:, 1)>bin_center(i), 1);
        speed_nearby = bin_speed([last_si, next_si], :);
        sp_lm = fitlm(speed_nearby(:, 1), speed_nearby(:, 2));
        center_sp = predict(sp_lm, bin_center(i));
        Speed(i) = center_sp;
    
    end

    Spike_Data_Bin{k} = Spikes;
    Position_Data_Bin{k} = Position;
    Speed_Data_Bin{k} = Speed;

end

%% Save the data
folder = 'D:\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
save(fullfile(folder, "Spike_Data_Binned.mat"), "Spike_Data_Bin", "-v7");
save(fullfile(folder, "Speed_Data_Binned.mat"), "Position_Data_Bin", "-v7");
save(fullfile(folder, "Position_Data_Binned.mat"), "Speed_Data_Bin", "-v7");


