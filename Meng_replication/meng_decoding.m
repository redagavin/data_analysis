% This is code reproduction of Meng paper using 04_tmaze1


%% Extract behavioral data from raw data, clean up
clc;
clear all;
folder = 'D:\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% parameter
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


nSD_frmap_smth = 3;
bin_sz_cm = 2;
dsp_delay = 984;

s_fname_bhv = 'bhv2.mat';
s_fname_suff_plf = '_plf.mat';
s_fname_csv = 'Shuo_linear_track_plf.csv';

lcnt = 1;
list_out = cell(1,1);

%% Load and pre-process behavioral data

nvt_file = fullfile(folder, PARAM.nvt_fname_in);

[TSusec, Xpos1, Ypos1, Ang1] = NlxNvtLoadYflip(nvt_file, PARAM.pix2cm);
[~, Xpos1, Ypos1] = BehavCleanUp(TSusec, Xpos1, Ypos1, PARAM);
BEHAV = NEW_BehavPreProc(TSusec, Xpos1, Ypos1, Ang1, PARAM.nSDsm_pos, PARAM.nSDsm_vel);

%% Linearize behavioral data

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

zz1= [];
tz1 = [];
vz1 = [];
tt = BEHAV.pos_ts_usec';
vv = BEHAV.vel_cmsec';

BEHAV.pos_ts_usec_z = tt;
BEHAV.pos_z_cm = zz;
BEHAV.vel_cmsec_z = vv;

BEHAV.outline_Z = [min(BEHAV.pos_z_cm) max(BEHAV.pos_z_cm)];
BEHAV.outline_legend_Z{1} = 'Zmin';
BEHAV.outline_legend_Z{2} = 'Zmax';

BEHAV.pos_x_cm = BEHAV.pos_z_cm';
BEHAV.pos_y_cm = ones(size(BEHAV.pos_x_cm));
BEHAV.outline = [BEHAV.outline_Z, 1, 1];

BEHAV.vel_cmsec = BEHAV.vel_cmsec_z';
BEHAV.pos_ts_usec = BEHAV.pos_ts_usec_z';

PARAM.bhv_cx = BEHAV.outline_Z(1) + (BEHAV.outline_Z(2) - BEHAV.outline_Z(1))/2;
PARAM.bhv_cy = BEHAV.outline(3) + (BEHAV.outline(4) - BEHAV.outline(3))/2;
PARAM.bhv_hsz_x = (BEHAV.outline_Z(2) - BEHAV.outline_Z(1))/2;
PARAM.bhv_hsz_y = (BEHAV.outline(4) - BEHAV.outline(3))/2;


PARAM.bhv_fname_out = fullfile(folder, 'bhv_linearized.mat');
save(PARAM.bhv_fname_out, 'PARAM', 'BEHAV');

%% Generate Position_Data_Processed for decoding

Position_Data(:,1)=BEHAV.pos_ts_usec_z/1e6;
Position_Data(:,2)=BEHAV.pos_x_cm;
Position_Data(:,3)=BEHAV.pos_y_cm;
Position_Data(:,4)=0;
TMP_BEHAV = NEW_BehavPreProc(Position_Data(:,1)*1e6, Position_Data(:,2), Position_Data(:,3), BEHAV.head_angles, PARAM.nSDsm_pos, PARAM.nSDsm_vel);
Position_Data(:,5)=TMP_BEHAV.vel_cmsec;
Position_Data(:,6)=0;
Position_Data(2:end,7)=diff(Position_Data(:,1));
Position_Data(1,7)=Position_Data(2,7);
Position_Data(Position_Data(:,7)>60,7)=0;
save(fullfile(folder, 'Position_Data_Processed.mat'),'Position_Data');

%% Extract Spike Data
clc;
clear all;
folder = 'D:\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% parameter

s_fname_csv = 'place_field_linearized.csv';%要填写的csv文件的文件名
s_fname_suff_plf = '_place_field_linearized.mat';%要描述place field的文件的文件名的后半部分
s_fname_bhv = 'bhv_linearized.mat';%描述behavior的文件的文件名
dsp_delay = 984;%系统产生的延迟常数
bin_sz_cm = 2;%每一帧的距离

PARAM.theta_band(1)=6;
PARAM.theta_band(2)=10;
PARAM.decim_r = 8;

lcnt = 1;
list_out = cell(1,1);

%% load and pre-process

BHV = load(fullfile(folder,s_fname_bhv));
TS_TRIAL = [BHV.BEHAV.pos_ts_usec(1) BHV.BEHAV.pos_ts_usec(end)];

Spike_Data = [];
current_cell_id = 0;   
Excitatory_Neurons = [];
Inhibitory_Neurons = [];
Tetrode_Cell_IDs = [];

ntt_file_list = LoadTargetFlist(folder,'cut_TT*.ntt');

for fid = 1:numel(ntt_file_list)%对每个ntt文件进行处理，下面的cell_ts, cell_ids都是该ntt文件中的，所有的参数，比如waveform, spike train, place fields, cell type，这些函数都是处理一个ntt文件的，或者说处理包含在该ntt文件内的所有细胞
  fprintf('\tProcess file: %s\n', ntt_file_list{fid});
  
  [s, s_ntt] = fileparts(ntt_file_list{fid});
  tt_num = s_ntt(isstrprop(s_ntt,'digit'));
  NCS_file_name = LoadTargetFlist(folder, ['*', num2str(tt_num),'_FSI*.ncs']);   
  LFP_WBAND = [];
  LFP_TS_USEC = [];
  LFP_srate = [];     
  LFP_theta = [];      
  [LFP_WBAND, LFP_TS_USEC] = NlxNcsGetAll(NCS_file_name{1});
  if NlxNcsIsInverted(NCS_file_name{1})
     LFP_WBAND = -LFP_WBAND;
  end
  LFP_srate = NlxNcsGetSFreq(NCS_file_name{1});
  LFP_theta = eegfilt(LFP_WBAND, LFP_srate, PARAM.theta_band(1),PARAM.theta_band(2) );

  
  [all_ts, all_cids, all_wfs, all_fets] = NlxGetSpikesAll( ntt_file_list{fid} );
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
  

  save( strcat(ntt_file_list{fid}(1:end-4),s_fname_suff_plf), ...
     'cell_ts', 'cell_ids', 'WFP', 'STP', ...
     'PLF', 's_*', 'CT', 'bin_sz_cm', 'LFP_TS_USEC', 'LFP_WBAND', 'LFP_theta', 'LFP_srate' ...
  );

  % Spike_Data_Processed
  for i = 1:size(cell_ids, 1)
      current_cell_id = current_cell_id + 1;
      current_spikes = cell_ts{i};
      current_cell_ids = current_cell_id * ones(numel(current_spikes), 1);
      current_spike_data = [current_spikes, current_cell_ids];
      Spike_Data = [Spike_Data; current_spike_data];

      
      if strcmp(CT{i}, 'INN_1')
          Inhibitory_Neurons = [Inhibitory_Neurons; current_cell_id];
      else
          Excitatory_Neurons = [Excitatory_Neurons; current_cell_id];
      end
      Tetrode_Cell_IDs = [Tetrode_Cell_IDs;[current_cell_id,str2num(tt_num), cell_ids{i}] ];
  end
end

[sort_1, sort_2] = sort(Spike_Data(:,1));
Spike_Data(:,1) = sort_1/1e6;
temp_Spike_Data2 = Spike_Data(:,2);
Spike_Data(:,2) = temp_Spike_Data2(sort_2);
save(fullfile(folder, 'Spike_Data_Processed.mat'),'Spike_Data', 'Excitatory_Neurons', 'Inhibitory_Neurons', 'Tetrode_Cell_IDs');

%% Spike and Position data integration
clc;
clear all;
folder = 'D:\memolab\data_analysis\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% load position and spike data

load Position_Data_Processed
load Spike_Data_Processed

%% Integration

Minimum_Time_Difference=0.1;

Matched_Spike_Data=zeros(size(Spike_Data,1),5);
parfor N=1:size(Spike_Data,1)
    Index=abs(Position_Data(:,1)-Spike_Data(N,1))==min(abs(Position_Data(:,1)-Spike_Data(N,1)));
    if sum(Index)>1
        Index=find(abs(Position_Data(:,1)-Spike_Data(N,1))==min(abs(Position_Data(:,1)-Spike_Data(N,1))),1,'first');
    end
    Matched_Spike_Data(N,:)=Position_Data(Index,1:5);
end
clear Index;

Spike_Information=Matched_Spike_Data(:,2:5);
Time_Difference=abs(Spike_Data(:,1)-Matched_Spike_Data(:,1));
Index=find(Time_Difference<=Minimum_Time_Difference);
fprintf('%d out of %d spikes were within %d seconds of the closest position timepoint. Max time difference was %d\n',length(Index),size(Spike_Data,1),Minimum_Time_Difference,max(Time_Difference)),
Spike_Data=Spike_Data(Time_Difference<=Minimum_Time_Difference,:);
Spike_Information=Spike_Information(Time_Difference<=Minimum_Time_Difference,:);
save(fullfile(folder, 'Spike_Information.mat'),'Spike_Information');