clc;
clear all;

%% loading and parameter
load Spike_Data_Separated.mat
load Position_Data_Separated.mat

bin_size = 0.05;
stride = 0.05;

Gausswin_size = 1;
Gausswin_sd = 0.05;

%% Bin
Spike_Data_Bin = {};
Position_Data_Bin = {};
Speed_Data_Bin = {};

Spike_Data_Bin_Right = {};
Position_Data_Bin_Right = {};
Speed_Data_Bin_Right = {};

Spike_Data_Bin_Left = {};
Position_Data_Bin_Left = {};
Speed_Data_Bin_Left = {};


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

    for i = 1:size(Spikes, 2)

        %kernel smoothing
        gauss_window = Gausswin_size/bin_size;
        gauss_SD = Gausswin_sd/bin_size;
        gk = gausskernel(gauss_window, gauss_SD)';

        Spikes(:,i) = conv2(Spikes(:,i),gk,'same');
        
    end

    % make spike count into firing rate
    Spikes = Spikes/bin_size;

    % normalize
    Spikes = zscore(Spikes, 0, 1);

    %pca
    [coeff, score, latent, tsquared, explained] = pca(Spikes);
    Spikes = Spikes * coeff(:, 1:40);

    Spike_Data_Bin{k} = Spikes;
    Position_Data_Bin{k} = Position;
    Speed_Data_Bin{k} = Speed;

end

[Spike_Data_Bin, trans] = hyperalign(Spike_Data_Bin{:});

% separate them into left and right groups

criterion = cellfun(@(x) any(x > 250), Position_Data_Bin);
Position_Data_Bin_Left = Position_Data_Bin(criterion);
Speed_Data_Bin_Left = Speed_Data_Bin(criterion);
Spike_Data_Bin_Left = Spike_Data_Bin(criterion);
Position_Data_Bin_Right = Position_Data_Bin(~criterion);
Speed_Data_Bin_Right = Speed_Data_Bin(~criterion);
Spike_Data_Bin_Right = Spike_Data_Bin(~criterion);


%% Save the data
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
save(fullfile(folder, "Spike_Data_Binned_Left.mat"), "Spike_Data_Bin_Left", "-v7");
save(fullfile(folder, "Speed_Data_Binned_Left.mat"), "Speed_Data_Bin_Left", "-v7");
save(fullfile(folder, "Position_Data_Binned_Left.mat"), "Position_Data_Bin_Left", "-v7");
save(fullfile(folder, "Spike_Data_Binned_Right.mat"), "Spike_Data_Bin_Right", "-v7");
save(fullfile(folder, "Speed_Data_Binned_Right.mat"), "Speed_Data_Bin_Right", "-v7");
save(fullfile(folder, "Position_Data_Binned_Right.mat"), "Position_Data_Bin_Right", "-v7");


