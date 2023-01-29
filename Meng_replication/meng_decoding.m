% This is code reproduction of Meng paper using 04_tmaze1


%% Extract behavioral data from raw data, clean up
clc;
clear all;
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

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
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

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
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

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

%% Place Field
clc;
clear all;
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% parameter and loading

load Position_Data_Processed
load Spike_Data_Processed
load Spike_Information

Bin_Size = 1;Velocity_Cutoff = 10;Firing_Rate_Cutoff = 0;Analyze_Linear=1;

%% Calculate Place Field

Pos_x_value = Position_Data(:,2);
Spike_x_value = Spike_Information(:,1);

Spike_Information(:,1)=Spike_Information(:,1)-min(Position_Data(:,2))+0.001;
Spike_Information(:,2)=Spike_Information(:,2)-min(Position_Data(:,3))+0.001;
Position_Data(:,2)=Position_Data(:,2)-min(Position_Data(:,2))+0.001;
Position_Data(:,3)=Position_Data(:,3)-min(Position_Data(:,3))+0.001;

Time_Change=diff(Position_Data(:,1));
Time_Change(end+1)=Time_Change(end);
Time_Change(Time_Change<=0)=min(Time_Change(Time_Change>0))/10;

X_Movement=diff(Position_Data(:,2));
X_Movement(end+1)=X_Movement(end);
Y_Movement=diff(Position_Data(:,3));
Y_Movement(end+1)=Y_Movement(end);
Position_Change=[X_Movement,Y_Movement];

Velocity=Position_Data(:,5);

Spike_Movement=zeros(size(Spike_Data,1),2);
parfor N=1:size(Spike_Data,1)
    Index=abs(Position_Data(:,1)-Spike_Data(N,1))==min(abs(Position_Data(:,1)-Spike_Data(N,1)));
    if sum(Index)>1
        Index=find(abs(Position_Data(:,1)-Spike_Data(N,1))==min(abs(Position_Data(:,1)-Spike_Data(N,1))),1,'first');
    end
    Spike_Movement(N,:)=[X_Movement(Index),Y_Movement(Index)];
end

xposition_Cutoff=34;
Index=Velocity>=Velocity_Cutoff & Pos_x_value>=xposition_Cutoff;
Position_Data=Position_Data(Index,:);
Time_Change=Time_Change(Index,:);
Position_Change=Position_Change(Index,:);
Position_Data(:,2)=ceil(Position_Data(:,2)/Bin_Size);
Position_Data(:,3)=ceil(Position_Data(:,3)/Bin_Size);
Index=Spike_Information(:,4)>=Velocity_Cutoff & Spike_x_value>=xposition_Cutoff;
Spike_Data=Spike_Data(Index,:);
Spike_Information=Spike_Information(Index,:);
Spike_Information(:,1)=ceil(Spike_Information(:,1)/Bin_Size);
Spike_Information(:,2)=ceil(Spike_Information(:,2)/Bin_Size);
Spike_Movement=Spike_Movement(Index,:);
clear Index;

Time_In_Position=zeros(max([max(Position_Data(:,3)),max(Spike_Information(:,2))]),max([max(Position_Data(:,2)),max(Spike_Information(:,1))]));
for N=1:size(Position_Data,1)
    Time_In_Position(Position_Data(N,3),Position_Data(N,2))=Time_In_Position(Position_Data(N,3),Position_Data(N,2))+Time_Change(N,1);
end

Spikes_In_Position=zeros(max([max(Position_Data(:,3)),max(Spike_Information(:,2))]),max([max(Position_Data(:,2)),max(Spike_Information(:,1))]),max(Spike_Data(:,2)));
for N=1:size(Spike_Data,1)
    Spikes_In_Position(Spike_Information(N,2),Spike_Information(N,1),Spike_Data(N,2))=Spikes_In_Position(Spike_Information(N,2),Spike_Information(N,1),Spike_Data(N,2))+1;
end

if Analyze_Linear
    if (max(Position_Data(:,2))-min(Position_Data(:,2)))>=(max(Position_Data(:,3))-min(Position_Data(:,3))) %if the track is horizontal
        Linear_Spikes_In_Position=sum(Spikes_In_Position,1);%如果track是横着的话，那么把所有y方向的spikes合并
        Linear_Spikes_In_Position=permute(Linear_Spikes_In_Position,[2,3,1]);%本来是%x,y,细胞index；现在变成y,细胞index,x
        Linear_Time_In_Position=sum(Time_In_Position,1);%如果track是横着的话，那么把所有y方向的time合并
        Linear_Time_In_Position=permute(Linear_Time_In_Position,[2,1]);
    elseif (max(Position_Data(:,2))-min(Position_Data(:,2)))<(max(Position_Data(:,3))-min(Position_Data(:,3))) %if the track is vertical
        Linear_Spikes_In_Position=sum(Spikes_In_Position,2);
        Linear_Spikes_In_Position=permute(Linear_Spikes_In_Position,[1,3,2]);%本来是%x,y,细胞index；现在变成x,细胞index,y
        Linear_Time_In_Position=sum(Time_In_Position,2);
    end
    Linear_Rate_In_Position=zeros(size(Linear_Spikes_In_Position,1),size(Linear_Spikes_In_Position,2));
    for N=1:size(Linear_Spikes_In_Position,2)
        Linear_Rate_In_Position(:,N)=Linear_Spikes_In_Position(:,N)./Linear_Time_In_Position;
    end
    Linear_Rate_In_Position(isnan(Linear_Rate_In_Position))=0;
    Linear_Rate_In_Position(isinf(Linear_Rate_In_Position))=0;
    Smoothed_Linear_Firing_Rate=Linear_Rate_In_Position;
    Filter=fspecial('gaussian',[20 1],2); %This filter smooths the firing rates with a gaussian kernel with a St.Dev. of 2 bins (4cm), out to 10 bins on either side (20cm)
    for N=1:size(Linear_Rate_In_Position,2)
        Smoothed_Linear_Firing_Rate(:,N)=filtfilt(Filter,1,Linear_Rate_In_Position(:,N));
    end
    Smoothed_Linear_Firing_Rate(isnan(Smoothed_Linear_Firing_Rate))=0;
    Smoothed_Linear_Firing_Rate(Smoothed_Linear_Firing_Rate<0)=0;
    for Z=1:size(Smoothed_Linear_Firing_Rate,2)
        if max(max(max(Smoothed_Linear_Firing_Rate(:,Z))))<Firing_Rate_Cutoff
            Smoothed_Linear_Firing_Rate(:,Z)=0;
        end
    end
    Field_Data_Linear=Smoothed_Linear_Firing_Rate;

    Firing_Rate_Peaks=zeros(max(Spike_Data(:,2)),4);
    Firing_Rate_Peaks(:,1)=1:max(Spike_Data(:,2));
    for N=1:size(Firing_Rate_Peaks)
        Firing_Rate_Peaks(N,2)=find(Field_Data_Linear(:,N)==max(Field_Data_Linear(:,N)),1,'first');
    end
    
    save(fullfile(folder, 'Field_Data.mat'),'Field_Data_Linear','Firing_Rate_Peaks');
end

%% Global Theta Filtered LFP
clc;
clear all;
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% parameter and loading

s_target_site = 'CA1';

P = struct();
P.decim_r = 8;
P.trg_band = [4 12];

BHV = load(fullfile(folder, 'bhv_linearized.mat'));
NCS_file_list = LoadTargetFlist(folder, '*_FSI*.ncs');
LFP_PWR = zeros(size(NCS_file_list));

S_BAND = [4  12];
N_BAND = [45 55];
snr_thr_dB = 6;

%% LFP

for fid = 1:numel(NCS_file_list)
    [LFP_WBAND, LFP_TS_USEC] = NlxNcsGetAll(NCS_file_list{fid});
    if NlxNcsIsInverted(NCS_file_list{fid})
        LFP_WBAND = -LFP_WBAND;
    end
    lfp_srate = NlxNcsGetSFreq(NCS_file_list{fid});

    LFP_WBAND = decimate(LFP_WBAND, P.decim_r);
    lfp_srate = lfp_srate / P.decim_r;

    LFP_FLT = eegfilt(LFP_WBAND, lfp_srate, P.trg_band(1), P.trg_band(2));
    LFP_PWR(fid) = sum(abs(hilbert(LFP_FLT)));

    [ PERC_SAT(fid), PERC_BAD_SNR(fid) ] = NlxNcsQuality( NCS_file_list{fid}, S_BAND, N_BAND, snr_thr_dB );
end

[~, best_ncs_idx] = min(PERC_BAD_SNR);

LFP_Filename = NCS_file_list{best_ncs_idx};
LFP_Frequency=Nlx2MatCSC(LFP_Filename,[0 0 1 0 0],0,3,1);
LFP_Header=Nlx2MatCSC(LFP_Filename,[0 0 0 0 0],1,1,0);
for Header_Line=1:length(LFP_Header)
    Header_Info=cell2mat(LFP_Header(Header_Line));
    if length(Header_Info)>12
        if strcmp(Header_Info(2:12),'-ADMaxValue')
            Max_Value=str2double(Header_Info(13:end-1));
        end
        if strcmp(Header_Info(2:12),'-InputRange')
            Max_Range=str2double(Header_Info(13:end));
        end
    end
end

[LFP_Times,LFP_Samples]=Nlx2MatCSC(LFP_Filename,[1 0 0 0 1],0,1);
LFP_Times=LFP_Times/1000000; 
LFP_Samples=LFP_Samples(:)*(Max_Range/Max_Value); 
Times=zeros(512,size(LFP_Times,2));
for B=1:length(LFP_Times)-1
    Times(:,B)=LFP_Times(B)+((0:(511))*((LFP_Times(B+1)-LFP_Times(B))/(512)))';
end
Times(:,end)=LFP_Times(end)+((0:(511))*((LFP_Times(end)-LFP_Times(end-1))/(512)))';
clear B;
clear LFP_Times;
Times=Times(:);
LFP_Data=[Times,LFP_Samples];
LFP_Left=LFP_Data;    
      
s1=strfind(NCS_file_list{fid},'CSC');
s2 =strfind(NCS_file_list{fid},'_FSI');
LFP_Electrodes=str2num(NCS_file_list{1}(s1+3:s2-1));

save(fullfile(folder, 'LFP_Left.mat'),'LFP_Left','LFP_Frequency','LFP_Electrodes','-v7.3');

clear LFP_Left;
clear Times;
clear LFP_Samples;
clear Max_Range;
clear Max_Value;

Theta_Stop_Low=5;
Theta_Pass_Low=6;                % Some papers use as low as 4 Hz as the low cutoff
Theta_Pass_High=12;              % Some papers use as high as 12 Hz as the high cutoff
Theta_Stop_High=14;
Stop_Band_Attenuation_One=60;    % This was the default, I think.
Pass_Band=1;                     % This was the default, I think.
Stop_Band_Attenuation_Two=80;    % This was the default, I think.
Filter_Design_For_Theta=fdesign.bandpass(Theta_Stop_Low, Theta_Pass_Low, Theta_Pass_High, Theta_Stop_High, Stop_Band_Attenuation_One, Pass_Band, Stop_Band_Attenuation_Two, LFP_Frequency);
Theta_Filter=design(Filter_Design_For_Theta,'butter');  %'equiripple' and 'cheby1' and 'cheby2' also work, but 'butter' was generally faster and gave similar results

Theta_Filtered_LFP_Data=zeros(size(LFP_Data,1),4);
Theta_Filtered_LFP_Data(:,1)=LFP_Data(:,1);
Theta_Filtered_LFP_Data(:,2)=filter(Theta_Filter,LFP_Data(:,2));
clear LFP_Data;
Theta_Filtered_LFP_Data(:,2)=Theta_Filtered_LFP_Data(end:-1:1,2);   
Theta_Filtered_LFP_Data(:,2)=filter(Theta_Filter,Theta_Filtered_LFP_Data(:,2));
Theta_Filtered_LFP_Data(:,2)=Theta_Filtered_LFP_Data(end:-1:1,2);
Theta_Filtered_LFP_Data(:,3)=hilbert(Theta_Filtered_LFP_Data(:,2));
Theta_Filtered_LFP_Data(:,4)=(angle(Theta_Filtered_LFP_Data(:,3))*180/pi)+180;
Theta_Filtered_LFP_Data(:,3)=abs(Theta_Filtered_LFP_Data(:,3));
Theta_Gaussian_Filter=fspecial('gaussian',[round(3*(300/((1/LFP_Frequency)*1000))),1],round(300/((1/LFP_Frequency)*1000))); %Gaussian filter with a sigma of 300 ms
Theta_Filtered_LFP_Data(:,3)=filtfilt(Theta_Gaussian_Filter,1,Theta_Filtered_LFP_Data(:,3));
Theta_Filtered_LFP_Data(:,3)=zscore(Theta_Filtered_LFP_Data(:,3));
LFP_Left_Theta=Theta_Filtered_LFP_Data;

save(fullfile(folder,'LFP_Left_Theta.mat'),'LFP_Left_Theta','-v7.3');

%% Theta Oscillation Properties
clc;
clear all;
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% parameter and loading
load('Spike_Data_Processed.mat');
load( 'Position_Data_Processed.mat');
load( 'Field_Data.mat');
load('Spike_Information.mat');
load( 'LFP_Left_Theta.mat') ;
load('LFP_Left.mat');

Theta_Length_Min_Max=[70/1e3 250/1e3];
Decoding_Time_Window = 20/1e3;Decoding_Time_Advance = 5/1e3;
Theta_Sequence_Starting_Phase=190;

%% Calculate properties

Excitatory_Spike_Data=Spike_Data;
for N=1:length(Inhibitory_Neurons)
    Excitatory_Spike_Data=Excitatory_Spike_Data(Excitatory_Spike_Data(:,2)~=Inhibitory_Neurons(N),:);
end
clear N;
LFP_Theta = LFP_Left_Theta;
Theta=LFP_Theta(:,[1,4,3]);
Starting_Time=min(Spike_Data(:,1));
Ending_Time=max(Spike_Data(:,1));
Decoding_Times=(Starting_Time+(Decoding_Time_Window/2):Decoding_Time_Advance:Ending_Time-(Decoding_Time_Window/2))';

Decoding_Time_Info=zeros(size(Decoding_Times,1),14);
Decoding_Time_Info_With_Inhibitory_Neurons=zeros(size(Decoding_Times,1),14);

Max_Number_Of_Spikes_Per_Decoding_Window=200;
Decoding_Spike_Index=zeros(Max_Number_Of_Spikes_Per_Decoding_Window,length(Decoding_Times));

parfor N=1:size(Decoding_Times,1)
    Participating_Spikes=Excitatory_Spike_Data(Excitatory_Spike_Data(:,1)>=(Decoding_Times(N)-(Decoding_Time_Window/2)) & Excitatory_Spike_Data(:,1)<=(Decoding_Times(N)+(Decoding_Time_Window/2)),2);
    Participating_Spikes_With_Inhibitory_Neurons=Spike_Data(Spike_Data(:,1)>=(Decoding_Times(N)-(Decoding_Time_Window/2)) & Spike_Data(:,1)<=(Decoding_Times(N)+(Decoding_Time_Window/2)),2);
    Position_Information=Position_Data(find(abs(Position_Data(:,1)-Decoding_Times(N))==min(abs(Position_Data(:,1)-Decoding_Times(N))),1,'first'),2:6);
    Decoding_Time_Info(N,:)=[Position_Information,0,size(Participating_Spikes,1),length(unique(Participating_Spikes)),0,0,0,0,0,0];
    Decoding_Time_Info_With_Inhibitory_Neurons(N,:)=[Position_Information,0,size(Participating_Spikes_With_Inhibitory_Neurons,1),length(unique(Participating_Spikes_With_Inhibitory_Neurons)),0,0,0,0,0,0];
    if size(Participating_Spikes,1)<Max_Number_Of_Spikes_Per_Decoding_Window
        Participating_Spikes(end+1:Max_Number_Of_Spikes_Per_Decoding_Window,:)=0;
    elseif size(Participating_Spikes,1)>Max_Number_Of_Spikes_Per_Decoding_Window
        error('Too many spikes detected.  Change Max_Number_Of_Spikes_Per_Decoding_Window in IRFS_CALCULATE_THETA_OSCILLATION_PROPERTIES to a larger value!')
    end
    Decoding_Spike_Index(:,N)=Participating_Spikes;
end
clear Position_Information;
clear Participating_Spikes;
clear Participating_Spikes_With_Inhibitory_Neurons;
clear N;
clear Max_Number_Of_Spikes_Per_Decoding_Window

Current_Line=1;
while max(Decoding_Spike_Index(Current_Line,:))>0
    Current_Line=Current_Line+1;
end
Decoding_Spike_Index=Decoding_Spike_Index(1:(Current_Line-1),:);
clear Current_Line;

All_Data_To_Add=zeros(size(Decoding_Times,1),2);
parfor N=1:size(Decoding_Times,1)  
    Data_To_Add=Theta(find(abs(Theta(:,1)-Decoding_Times(N,1))==min(abs(Theta(:,1)-Decoding_Times(N,1))),1,'first'),2:3);
    All_Data_To_Add(N,:)=Data_To_Add;
end
Decoding_Time_Info(:,[6,11])=All_Data_To_Add;
Decoding_Time_Info_With_Inhibitory_Neurons(:,[6,11])=All_Data_To_Add;
clear N;
clear Data_To_Add;
clear All_Data_To_Add;

Current_Oscillation_ID_Number=0;
for N=2:size(Decoding_Time_Info,1)
    if Decoding_Time_Info(N-1,6)<Theta_Sequence_Starting_Phase && Decoding_Time_Info(N,6)>=Theta_Sequence_Starting_Phase 
        Current_Oscillation_ID_Number=Current_Oscillation_ID_Number+1;
    end
    Decoding_Time_Info(N,14)=Current_Oscillation_ID_Number;
end
clear N;
clear Current_Oscillation_ID_Number;

Velocity_Cutoff = 5;xposition_Cutoff=34; Initial_Variables.Decoding_Time_Advance=5/1e3;

for N=1:max(Decoding_Time_Info(:,14))
    Index=find(Decoding_Time_Info(:,14)==N);
    if ~isempty(Index)
        Current_Decoding_Time_Info=Decoding_Time_Info(Index,:);
        if isempty(find(Current_Decoding_Time_Info(:,4)<Velocity_Cutoff,1)) && isempty(find(Current_Decoding_Time_Info(:,1)<xposition_Cutoff, 1))
            Decoding_Time_Info(Index,9)=1;
        end
        Oscillation_Duration=size(Current_Decoding_Time_Info,1)*Initial_Variables.Decoding_Time_Advance;
        Current_Oscillation_Theta_Phases=Current_Decoding_Time_Info(:,6);
        Current_Oscillation_Theta_Phases(Current_Oscillation_Theta_Phases<=Theta_Sequence_Starting_Phase)=Current_Oscillation_Theta_Phases(Current_Oscillation_Theta_Phases<=Theta_Sequence_Starting_Phase)+360;
        if isempty(find(diff(Current_Oscillation_Theta_Phases)<0,1)) && Oscillation_Duration>=Theta_Length_Min_Max(1) && Oscillation_Duration<=Theta_Length_Min_Max(2)
            Decoding_Time_Info(Index,10)=1;
        end
    end
    clear Index;
    clear Oscillation_Duration;
    clear Current_Oscillation_Theta_Phases;
    clear Current_Decoding_Time_Info;
    clear Unimodal_Window_Spike_Count;
    clear Bimodal_Window_Spike_Count;
end
clear N;

Decoding_Window_Index=find(Decoding_Time_Info(:,7)>0 & Decoding_Time_Info(:,9)==1 & Decoding_Time_Info(:,10)==1);
Decoding_Window_Index_With_Inhibitory_Neurons=find(Decoding_Time_Info_With_Inhibitory_Neurons(:,7)>0 & Decoding_Time_Info_With_Inhibitory_Neurons(:,9)==1 & Decoding_Time_Info_With_Inhibitory_Neurons(:,10)==1);

save('Decoding_Time_Info','Decoding_Times','Decoding_Time_Info','Decoding_Time_Info_With_Inhibitory_Neurons','Starting_Time','Ending_Time','Decoding_Window_Index','Decoding_Window_Index_With_Inhibitory_Neurons','-v7.3');
save('Decoding_Spike_Index_All_Decoding_Windows','Decoding_Spike_Index','-v7.3')

Decoding_Spike_Index=Decoding_Spike_Index(:,Decoding_Window_Index);
Position_Data_For_Shuffles=Decoding_Time_Info(Decoding_Window_Index,[1,2,5]);
save('Decoding_Spike_Index','Decoding_Spike_Index','Position_Data_For_Shuffles','-v7.3')

%% Decode linear theta sequences
clc;
clear all;
folder = '..\Tmaze_Data\252-1375\2018-01-07_15-14-54\04_tmaze1';

%% Parameter and Loading

load Decoding_Time_Info.mat
load Decoding_Spike_Index.mat
load Field_Data_Light.mat
load Position_Data_Processed
load Field_Data

Decoding_Time_Window = 20/1e3;Decoding_Time_Advance = 5/1e3;

%% decoding

Field_Data_Linear2=Field_Data_Linear;

for N=1:size(Field_Data_Linear2,2)
    Field=Field_Data_Linear2(:,N);
    Minimum=min(Field(Field>0));

    if ~isempty(Minimum)
        Field(Field<=0)=Minimum/10;
    else
        Field(:)=1;
    end

    Field_Data_Linear2(:,N)=Field;

    clear Minimum;
    clear Minimum_In;
    clear Minimum_Out;
    clear Field;
    clear Field_In;
    clear Field_Out;
end
clear N;

if max(Position_Data(:,2))-min(Position_Data(:,2))>max(Position_Data(:,3))-min(Position_Data(:,3))
    Position_Column=1;
else
    Position_Column=2;
end
Movement_Direction=diff(Decoding_Time_Info(:,Position_Column));
Movement_Direction(end+1,1)=Movement_Direction(end,1);
if Movement_Direction(1,1)==0
    Movement_Direction(1,1)=Movement_Direction(find(Movement_Direction>0,1,'first'));
end
for N=2:length(Movement_Direction)
    if Movement_Direction(N,1)==0
        Movement_Direction(N,1)=Movement_Direction(N-1,1);
    end
end
Movement_Direction(Movement_Direction>=0)=1;
Movement_Direction(Movement_Direction<0)=-1;

Mid_Point=round(size(Field_Data_Linear2,1)/2);

Decoded_Data_Linear_Raw=zeros(size(Field_Data_Linear2,1),length(Decoding_Window_Index),3); 
Decoded_Data_Linear=zeros(size(Field_Data_Linear2,1),length(Decoding_Window_Index));

for Current_Decoding_Window=1:length(Decoding_Window_Index)
    
    Subset_Spike_Data=Decoding_Spike_Index(:,Current_Decoding_Window);
    Subset_Spike_Data=Subset_Spike_Data(Subset_Spike_Data>0);
        
    Decoded_Matrix=prod(Field_Data_Linear2(:,Subset_Spike_Data),2).*exp(-Decoding_Time_Window*sum(Field_Data_Linear2,2));
    if isinf(max(max(Decoded_Matrix))) 
        Divider=1;
        while isinf(max(max(Decoded_Matrix)))
            Decoded_Matrix=prod((Field_Data_Linear2(:,Subset_Spike_Data)/(2^Divider)),2).*exp(-Decoding_Time_Window*sum((Field_Data_Linear2/(2^Divider)),2));
            Divider=Divider+1;
        end
    end
    Decoded_Matrix(Decoded_Matrix<0)=0;
    if max(max(Decoded_Matrix))>0
        Decoded_Matrix=Decoded_Matrix/sum(sum(Decoded_Matrix));
    end

    Bin_Size = 1;

    Location=min([size(Field_Data_Linear2,1),max([1,round(Decoding_Time_Info(Decoding_Window_Index(Current_Decoding_Window),Position_Column)/Bin_Size)])]);
    Direction=Movement_Direction(Decoding_Window_Index(Current_Decoding_Window),1);
    Translated_Image=imtranslate(Decoded_Matrix,[0,Mid_Point-Location],'FillValues',0);

    Decoded_Data_Linear(:,Current_Decoding_Window)=Translated_Image;
end

clear Current_Decoding_Window;
clear Subset_Spike_Data;
clear Divider;
clear Location;
clear Direction;
clear Translated_Image;
clear Decoded_Matrix;
clear Mid_Point;
clear Field_Data_Linear2;
clear Position_Column;
clear Decoded_Matrix;
clear Sum;
save(fullfile(folder, 'Decoded_Linear_Data_And_Sequences.mat'),'Movement_Direction','Decoded_Data_Lin*','Decoding_Time_Window','Decoding_Time_Advance','Decoding_Window_Index','-v7.3');
clear Movement_Direction;

