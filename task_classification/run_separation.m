clc;
clear all;

%% loading and parameter
load Spike_Data_Processed.mat
load Position_Data_Processed.mat

%% separation
Spike_Data_Sep = {};
Position_Data_Sep = {};

page = 1;
start = 1;
finish = 1;
for i = 1:length(Position_Data)
    if i == length(Position_Data)
        continue;
    elseif Position_Data(i+1,2) > 35 && Position_Data(i,2) <= 35
        start = i;
    elseif Position_Data(i, 2) > 105 && Position_Data(i+1, 2) < 35
        finish = i;
        if length(Position_Data(start:finish, :)) <= 300
            continue;
        end
        Position_Data_Sep{page} = Position_Data(start:finish, :);
        start_time = Position_Data(start, 1);
        finish_time = Position_Data(finish, 1);
        Spike_Data_Sep{page} = Spike_Data(Spike_Data(:,1)>=start_time & Spike_Data(:, 1)<=finish_time, :);
        page = page + 1;
    end
end

%% Save the data
folder = 'D:\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
save(fullfile(folder, "Spike_Data_Separated.mat"), "Spike_Data_Sep");
save(fullfile(folder, "Position_Data_Separated.mat"), "Position_Data_Sep")
