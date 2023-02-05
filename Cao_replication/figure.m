clc;
clear all;

%% load data
load('bin_data.mat');


vidx = find(Speed>0.01);    % speed threshold 5 cm/s
rpos_init = normalize(Position(vidx,:),1,'range');
rpos_init = rpos_init(9643:9920,:);

rspk = Spikes(vidx,:);
rspk = rspk(9643:9920,:);

%% decoding
deme = 'vb';
mnames={'Spikes', 'Spikes (+history)'};

median_error = zeros(0);
DecodedPos = zeros(0);

trackLen = 290;
OB=cell(0); plotline = 2;
rpos = rpos_init(2:end);
OB{1} = zscore(rspk(2:end,:));

OB{2} = zscore([rspk(2:end,:),rspk(1:end-1,:)]);


% OLE
runlist = 1:length(OB);
for im=runlist
    obs = OB{im};
    % Generate basis...
    [basis,Xrbf] = get1Dbasis('vonmises',75,rpos*2*pi,100);

    pvec = linspace(0,2*pi,trackLen/2); % circular phase 
    [~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);

    % Cross-validated predictions
    mcvp.cv = getCVidx(size(obs,1),10); f=[];
    for it=1:mcvp.cv.nfoldcv
        if strcmp(deme,'ls')
            % least squared estimation  training
            w = obs(mcvp.cv.tr{it},:) \ Xrbf(mcvp.cv.tr{it},:);
        elseif strcmp(deme,'vb')
            % Variational Bayes training
            w = zeros(0);
            for j=1:size(Xrbf(mcvp.cv.tr{it},:),2)
                [w(:,j),V,~,~,~,~,~,~] = ...
                    VB_ARD_linear_regression(obs(mcvp.cv.tr{it},:), Xrbf(mcvp.cv.tr{it},j));
            end
        end
        % testing
        f(mcvp.cv.ts{it},:) = obs(mcvp.cv.ts{it},:) * (w*dbasis');
    end

    % Plotting and error calcs... 
    f = zscore(f,[],2);
    idx = 1:size(f,1);
    [~,maxpost] = max(f,[],2);
    DecodedPos(:,2*(im-1)+1) = trackLen*rpos(idx);
    DecodedPos(:,2*(im-1)+2) = trackLen*maxpost/length(pvec);

    % Circular error...
    errp = circ_dist(     rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
    errm = circ_dist(2*pi-rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
    % Save results
    errpAll(:,im) = errp;
    errmAll(:,im) = errm;
    % non-directional error
    [err_nd_All(:,im),mini] = (min([abs(errmAll(:,im)) abs(errpAll(:,im))]'));
    pos_hat(:,im) = maxpost'/length(pvec)*2*pi;
end

% Median error in cm with s.e.
for im=runlist
    bootstat = bootstrp(500,'median',abs(err_nd_All(:,im)));
    median_error(im,:) = [mean(bootstat),std(bootstat)];
end

rpos = rpos*trackLen;
pos_hat = round(pos_hat*trackLen/(2*pi));
save('Figure1bc_results_data','trackLen','mnames','rpos','pos_hat','err_nd_All','median_error');

%% Plot
load('Figure1bc_results_data');
plotwindow = [-0.1,30];
timestamps = [1:length(rpos)]*0.1;
idx = [0;diff(rpos)] > 1.5;
rpos(idx) = nan; 
ms  = 36;

figure('Color','white');
plot(timestamps,rpos,'k','LineWidth',3); hold on;
s = scatter(timestamps,pos_hat(:,1),ms,'r','filled');
s.MarkerFaceAlpha = 0.5;
xlim(plotwindow);
title("Clustered spikes")

