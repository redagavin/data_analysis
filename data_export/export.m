clc;
clear all;

%% loading nvt and export it as csv
nvt_file = 'C:\Users\yzh_x\Desktop\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1\VT1_fixed.nvt';
pix2cm = 5.4135;

[TSusec, Xpos1, Ypos1, Ang1] = NlxNvtLoadYflip(nvt_file, pix2cm);
m_nvt = [TSusec, Xpos1, Ypos1, Ang1];
writematrix(m_nvt, 'nvt.csv');

%% loading ntt and export it as csv
ntt_folder = 'C:\Users\yzh_x\Desktop\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
ntt_file_list = LoadTargetFlist(ntt_folder, 'cut_TT*.ntt');

for f = 1:numel(ntt_file_list)
    [all_ts, all_cids, all_wfs, all_fets] = NlxGetSpikesAll(ntt_file_list{f});

    % all_fets' dimension is 8x(number of timestamps). transpose it.
    all_fets = transpose(all_fets);

    % all_wfs' dimension is 32x4x(number of timestamps). reshape and transpose it.
    wfs_size = size(all_wfs);
    all_wfs = reshape(all_wfs, [wfs_size(1)*wfs_size(2),wfs_size(3)]);
    all_wfs = transpose(all_wfs);

    m_ntt = [all_ts, all_cids, all_wfs, all_fets];
    [path, name, ext] = fileparts(ntt_file_list{f});
    writematrix(m_ntt, strcat(name, '.csv'));
end

%% loading ncs and export it as csv
ncs_folder = 'C:\Users\yzh_x\Desktop\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
ncs_file_list = LoadTargetFlist(ncs_folder, 'CSC*_FSI20.ncs');

for f = 1:numel(ncs_file_list)
    [LFP_WBAND, LFP_TS_USEC] = NlxNcsGetAll(ncs_file_list{f});
    m_ncs = [LFP_WBAND, LFP_TS_USEC];
    [path, name, ext] = fileparts(ncs_file_list{f});
    writematrix(m_ncs, strcat(name, '.csv'));
end