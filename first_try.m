clear all;
close all;
clc;
% load('/Users/bettylin2727/Downloads/AtlasViewer/Group/DataTree/AcquiredData/Snirf/Examples/neuro_run01.nirs', '-mat');
% load('/Users/bettylin2727/Downloads/AtlasViewer/Group/DataTree/AcquiredData/Snirf/Examples/neuro_run01.snirf');
[snirf_saved, snirf_loaded, nirs] = snirf_load_save('/Users/bettylin2727/Downloads/sub_liao.nirs');
nReps= 8;
time_perStim = 45;
sampling_freq = 50;
%% Remove bad channels according SCI
nb_channels = size(snirf_loaded.data.measurementList,2)/2;
sci = get_SCI(snirf_loaded.data, nb_channels);
disp (sci);
% disp(snirf_loaded.data.measurementList(1).sourceIndex);
% disp(snirf_loaded.data.measurementList(1).detectorIndex);

raw_intensity = snirf_loaded.data.dataTimeSeries;
time= snirf_loaded.data.time;
figure();
plot(time, raw_intensity(:,1), 'r'); hold on;
plot(time, raw_intensity(:,2), 'g'); hold on;
plot(time, raw_intensity(:,3), 'b');
xlabel('time (s)');
ylabel('raw intensity');


%%
dod = hmrR_Intensity2OD( snirf_loaded.data );
%dod = hmrR_Intensity2OD( snirf_loaded );
figure();
plot(time, dod.dataTimeSeries(:,1), 'r'); hold on;
plot(time, dod.dataTimeSeries(:,2), 'g'); hold on;
plot(time, dod.dataTimeSeries(:,3), 'b');
xlabel('time (s)');
ylabel('delta optical density');



%% Motion Artifact correction
% mlActMan, mlActAuto ??
% dod_MA_corrected = hmrR_MotionCorrectWavelet(dod, cell(length(dod),1),cell(length(dod),1) ,1.5);


% PCA remove one component: one source motion artifact
dod_MA_corrected = hmrR_MotionCorrectPCA(dod, [],[],[],[], 2);

figure();
plot(time, dod_MA_corrected.dataTimeSeries(:,1),'r'); hold on;
plot(time, dod_MA_corrected.dataTimeSeries(:,2),'g'); hold on;
plot(time, dod_MA_corrected.dataTimeSeries(:,3),'b'); 
xlabel('time (s)');
ylabel('delta optical density');
title('MA corrected');

%% Bandpass filter drift correction

data_bp = hmrR_BandpassFilt(dod_MA_corrected, 0.5, 0.01);
%% OD to Conc (apply MBLL)
dc = hmrR_OD2Conc( data_bp, snirf_loaded.probe, 1 );

%{
14 channels
3 concentration
8 reps
4 volumns
%}

%% Reduce extracerebral component with short channels
% ROI 1
for i = 0:2
    dc.dataTimeSeries(:, 1*3-i) = dc.dataTimeSeries(:, 1*3-i) - dc.dataTimeSeries(:, 13*3-i);
    dc.dataTimeSeries(:, 2*3-i) = dc.dataTimeSeries(:, 2*3-i) - dc.dataTimeSeries(:, 13*3-i);
    dc.dataTimeSeries(:, 5*3-i) = dc.dataTimeSeries(:, 5*3-i) - dc.dataTimeSeries(:, 13*3-i);

    % ROI3
    dc.dataTimeSeries(:, 3*3-i) = dc.dataTimeSeries(:, 3*3-i) - dc.dataTimeSeries(:, 13*3-i);
    dc.dataTimeSeries(:, 4*3-i) = dc.dataTimeSeries(:, 4*3-i) - dc.dataTimeSeries(:, 13*3-i);
    dc.dataTimeSeries(:, 8*3-i) = dc.dataTimeSeries(:, 8*3-i) - dc.dataTimeSeries(:, 13*3-i);

    % ROI5
    dc.dataTimeSeries(:, 14*3-i) = dc.dataTimeSeries(:, 14*3-i) - dc.dataTimeSeries(:, 13*3-i);
    dc.dataTimeSeries(:, 6*3-i) = dc.dataTimeSeries(:, 6*3-i) - dc.dataTimeSeries(:, 11*3-i);
    dc.dataTimeSeries(:, 9*3-i) = dc.dataTimeSeries(:, 9*3-i) - dc.dataTimeSeries(:, 11*3-i);

    % ROI7
    dc.dataTimeSeries(:, 7*3-i) = dc.dataTimeSeries(:, 7*3-i) - dc.dataTimeSeries(:, 11*3-i);
    dc.dataTimeSeries(:, 10*3-i) = dc.dataTimeSeries(:, 10*3-i) - dc.dataTimeSeries(:, 11*3-i);
    dc.dataTimeSeries(:, 12*3-i) = dc.dataTimeSeries(:, 12*3-i) - dc.dataTimeSeries(:, 11*3-i);

end
%% starting points

% starting time, unit is sec
start_idx = 1 + 50*  min([snirf_loaded.stim(1,1).data(1,1), snirf_loaded.stim(1,2).data(1,1), snirf_loaded.stim(1,3).data(1,1), snirf_loaded.stim(1,4).data(1,1)]);
% starting points for 4 stimuli
start_idx_A = 1.+ 50* snirf_loaded.stim(1,1).data(:,1);
start_idx_B = 1.+ 50* snirf_loaded.stim(1,2).data(:,1);
start_idx_C = 1.+ 50* snirf_loaded.stim(1,3).data(:,1);
start_idx_D = 1.+ 50* snirf_loaded.stim(1,4).data(:,1);

%% First set
A_1 = dc.dataTimeSeries(start_idx_A(1):start_idx_A(1)+ 45*50 -1 ,:);
B_1 = dc.dataTimeSeries(start_idx_B(1):start_idx_B(1)+ 45*50 -1 ,:);
C_1 = dc.dataTimeSeries(start_idx_C(1):start_idx_C(1)+ 45*50 -1 ,:);
D_1 = dc.dataTimeSeries(start_idx_D(1):start_idx_D(1)+ 45*50 -1 ,:);
%  Plot 14 channels
plot_14Ch(A_1, B_1, C_1, D_1);


%% Timeseries seperation
% devide vertically

cell_A = cell(1,8);
cell_B = cell(1,8);
cell_C = cell(1,8);
cell_D = cell(1,8);

sum_A = zeros(45*50, 42);
sum_B = zeros(45*50, 42);
sum_C = zeros(45*50, 42);
sum_D = zeros(45*50, 42);


for i = 1: nReps
    
    cell_A{1,i} = dc.dataTimeSeries(start_idx_A(i):start_idx_A(i)+ 45*50 -1 ,:);
    cell_B{1,i} = dc.dataTimeSeries(start_idx_B(i):start_idx_B(i)+ 45*50 -1 ,:);
    cell_C{1,i} = dc.dataTimeSeries(start_idx_C(i):start_idx_C(i)+ 45*50 -1 ,:);
    cell_D{1,i} = dc.dataTimeSeries(start_idx_D(i):start_idx_D(i)+ 45*50 -1 ,:);
    
    sum_A = sum_A + cell_A{1,i};
    sum_B = sum_B + cell_B{1,i};
    sum_C = sum_C + cell_C{1,i};
    sum_D = sum_D + cell_D{1,i};
    
end


avg_A = sum_A./nReps;
avg_B = sum_B./nReps;
avg_C = sum_C./nReps;
avg_D = sum_D./nReps;

%% Plot 14 channels
plot_14Ch(avg_A, avg_B, avg_C, avg_D);


%% 
printf('');