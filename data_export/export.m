clc;
clear all;

%% loading and preprocessing nvt and export it as csv

PARAM = struct();

PARAM.pix2cm = 5.4135;
PARAM.vel_thresh_cmsec = 2;
PARAM.lap_min_boundary = 25;
PARAM.lap_max_boundary = 110;


PARAM.nSDsm_pos = 0.01; 
PARAM.nSDsm_vel = 2.5;

PARAM.t_min_sec = 3;
PARAM.max_ret = 0.1;
PARAM.max_vel_cmsec = 40;
PARAM.nvt_fname_in = 'VT1_fixed.nvt';
PARAM.bhv_fname_out = ''; % to be filled in the loop below

PARAM.bhv_cx = 0;         % to be filled in the loop below
PARAM.bhv_cy = 0;         % to be filled in the loop below
PARAM.bhv_hsz_x = 0;     % half length of the track, cm % to be filled in the loop below
PARAM.bhv_hsz_y = 0;     % half width of the track, cm % to be filled in the loop below

PARAM.plf_num_colors = 64;
PARAM.plf_thr_type = 'perPFR';
PARAM.plf_thr_perPFR = 20;
PARAM.ds_thr = 10;
PARAM.sgolayfilt_k = 3;
PARAM.sgolayfilt_f = 31;
PARAM.interp1_method = 'nearest';

nvt_file = 'C:\Users\yzh_x\Desktop\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1\VT1_fixed.nvt';

[TSusec, Xpos1, Ypos1, Ang1] = NlxNvtLoadYflip(nvt_file, PARAM.pix2cm);
[~, Xpos1, Ypos1] = BehavCleanUp(TSusec, Xpos1, Ypos1, PARAM);
BEHAV = BehavPreProc(TSusec, Xpos1, Ypos1, Ang1, PARAM.nSDsm_pos, PARAM.nSDsm_vel);

xx = BEHAV.pos_x_cm;
yy = BEHAV.pos_y_cm;
for i = 1:numel(xx)
   if yy(i) <= 49 & yy(i) >= 43.5
       if xx(i)<111
           zz(i) = xx(i);
       else
           zz(i) = 111;
       end
   end
   if yy(i) > 49
       if xx(i)>= 108
           zz(i) = yy(i) - 49 + 111;
       elseif xx(i)< 108 & xx(i)>=34
           zz(i) = 111 + 81 - 49 + sqrt((81 - yy(i))^2 + (xx(i)-108)^2);
       else
           zz(i) = xx(i);
       end
   end
   if yy(i) < 43.5
       if xx(i)>= 108
           zz(i) = - yy(i) + 43.5 + 221.5;
       elseif xx(i)< 108 & xx(i)>=34
           zz(i) = 43.5 - 11.5 + sqrt((yy(i) - 11.5)^2 + (xx(i)-108)^2) + 221.5;
       else
           zz(i) = xx(i);
       end
   end
end

% m_nvt = [TSusec, Xpos1, Ypos1, Ang1];
% writematrix(m_nvt, 'nvt.csv');

%% loading ntt and export it as csv
ntt_folder = 'C:\Users\yzh_x\Desktop\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';
ntt_file_list = LoadTargetFlist(ntt_folder, 'cut_TT*.ntt');

dsp_delay = 984;

for f = 1:numel(ntt_file_list)
    [all_ts, all_cids, all_wfs, all_fets] = NlxGetSpikesAll(ntt_file_list{f});
    all_ts = all_ts - dsp_delay;

    cell_fet = SpkSelectFet(all_fets, all_cids, false);

    [cell_ts, cell_ids] = SpkSelectTs(all_ts, all_cids, false);

    cell_wfs = SpkSelectWf(all_wfs, all_cids, false);

    bad_cells = false(length(cell_ids),1);
    for ii = 1:length(cell_ids)
        if numel(cell_ts{ii}) < 50
            bad_cells(ii) = true;
        end
    end
    cell_fet(bad_cells) = [];
    cell_ids(bad_cells) = [];
    cell_ts(bad_cells)  = [];
    cell_wfs(bad_cells) = [];

    WFP = SpkWFormProp2(cell_ids, cell_wfs, 1, true);
    WFP_BEST_CH = [WFP(1).best_ch; vertcat(WFP.best_ch)];

    STP = SpkTrainProp2(cell_ids, cell_ts, cell_wfs, cell_fet, WFP_BEST_CH, TS_TRIAL, 1);

    cell_fet = cell_fet(2:end);
    cell_ids = cell_ids(2:end);
    cell_ts  = cell_ts(2:end);
    cell_wfs = cell_wfs(2:end);

    PLF = CalcPlaceFields(BHV.PARAM, BHV.BEHAV, cell_ts);

    CT  = CalcCellTypeWeak( WFP, STP );

    % all_fets' dimension is 8x(number of timestamps). transpose it.
%     all_fets = transpose(all_fets);

    % all_wfs' dimension is 32x4x(number of timestamps). reshape and transpose it.
%     wfs_size = size(all_wfs);
%     all_wfs = reshape(all_wfs, [wfs_size(1)*wfs_size(2),wfs_size(3)]);
%     all_wfs = transpose(all_wfs);

%     m_ntt = [all_ts, all_cids, all_wfs, all_fets];
%     [path, name, ext] = fileparts(ntt_file_list{f});
%     writematrix(m_ntt, strcat(name, '.csv'));
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