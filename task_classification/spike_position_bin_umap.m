clc;
clear all;

%% loading and parameter
load Spike_Data_Separated.mat
load Position_Data_Separated.mat

global Gausswin_size;
global Gausswin_sd;
Gausswin_size = 4;
Gausswin_sd = 1;

global bin_size;
bin_size = 1;
upper_bound = 219.5;
lower_bound = 35.5;

upper = lower_bound+1:bin_size:upper_bound;


%% Separate left and right
criterion = cellfun(@(x) any(x(:, 2) > 250), Position_Data_Sep);
Position_Data_Sep_Left = Position_Data_Sep(criterion);
Spike_Data_Sep_Left = Spike_Data_Sep(criterion);
Position_Data_Sep_Right = Position_Data_Sep(~criterion);
Spike_Data_Sep_Right = Spike_Data_Sep(~criterion);

%% process left position data

for i = 1:length(Position_Data_Sep_Left)

    Position_Data = Position_Data_Sep_Left{i};

    for j = 1:size(Position_Data, 1)
        if Position_Data(j, 2) >= 221
            Position_Data(j, 2) = Position_Data(j, 2) - 110.5;
        end
    end

    Position_Data_Sep_Left{i} = Position_Data;
end


%% Bin and UMAP

Spike_Data_Umap_Right = {};

Spike_Data_Umap_Left = {};

Spike_Data_Umap_Right = process_spike(Position_Data_Sep_Right, Spike_Data_Sep_Right, upper);

Spike_Data_Umap_Left = process_spike(Position_Data_Sep_Left, Spike_Data_Sep_Left, upper);


%% Save the data
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
save(fullfile(folder, "Spike_Data_Umap_Left.mat"), "Spike_Data_Umap_Left", "-v7");
save(fullfile(folder, "Spike_Data_Umap_Right.mat"), "Spike_Data_Umap_Right", "-v7");


%% function for Bin and UMAP
function umap_data = process_spike(Position_Data_Sep, Spike_Data_Sep, upper)
    
    global Gausswin_size;
    global Gausswin_sd;
    global bin_size;

    umap_data = {};
    
    % for each run-through, we bin the data
    for k = 1:length(Position_Data_Sep)
    
        % get the data
        Spike_Data = Spike_Data_Sep{k};
        Position_Data = Position_Data_Sep{k};
        
        % pre-fill the arrays with zeros
        % they will be filled with real data later in the loop
        Spikes = zeros(length(upper), max(unique(Spike_Data(:, 2))));
        
        for i = 1:length(upper)
        
            % get the upper bound and lower bound of the bin
            upper_b = upper(i);
            lower_b = upper_b - 1;
        
            % find the first position that is larger than the upper bound
            idx_l = find(Position_Data(:, 2)>upper_b, 1);
            
            % find the positions that are near the upper bound
            upper_nearby = Position_Data([idx_l-1, idx_l], [1, 2]);
    
            % find the time corresponding to the upper bound position
            upper_lm = fitlm(upper_nearby(:, 2), upper_nearby(:, 1));
            upper_time = predict(upper_lm, upper_b);
            
            % find the last position that is smaller than the lower bound and
            % the corresponding time step should be smaller than the upper
            % bound time step
            idx_s = find((Position_Data(:, 2)<lower_b) & ...
                (Position_Data(:, 1)<Position_Data(idx_l, 1)), 1, "last");
            
            % find the positions that are near the lower bound
            lower_nearby = Position_Data([idx_s, idx_s+1], [1, 2]);
    
            % find the time corresponding to the lower bound position
            lower_lm = fitlm(lower_nearby(:, 2), lower_nearby(:, 1));
            lower_time = predict(lower_lm, lower_b);
    
            % find all the spikes in the bin
            bin_spike = Spike_Data( ...
                (Spike_Data(:, 1)>=lower_time) & (Spike_Data(:, 1)<=upper_time), :);
            
            % process the spike data in the bin
            for j = bin_spike(:, 2)'
                Spikes(i, j) = Spikes(i, j) + 1;
            end
    
            % make spike count into firing rate
            Spikes(i, :) = Spikes(i, :)/(upper_time-lower_time);
    
        end
    
        for i = 1:size(Spikes, 2)
    
            %kernel smoothing
            gauss_window = Gausswin_size/bin_size;
            gauss_SD = Gausswin_sd/bin_size;
            gk = gausskernel(gauss_window, gauss_SD)';
    
            Spikes(:,i) = conv2(Spikes(:,i),gk,'same');
            
        end
    
        % normalize
        % Spikes = zscore(Spikes, 0, 1);
    
        % umap
        Spikes_umap = run_umap(Spikes,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3, 'verbose', 'text');
    
        umap_data{k} = Spikes_umap;
    
    end
end

