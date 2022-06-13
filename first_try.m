clear all;
close all;
clc;
% load('/Users/bettylin2727/Downloads/AtlasViewer/Group/DataTree/AcquiredData/Snirf/Examples/neuro_run01.nirs', '-mat');
% load('/Users/bettylin2727/Downloads/AtlasViewer/Group/DataTree/AcquiredData/Snirf/Examples/neuro_run01.snirf');
[snirf_saved, snirf_loaded, nirs] = snirf_load_save('/Users/bettylin2727/Downloads/AtlasViewer/Group/DataTree/AcquiredData/Snirf/Examples/neuro_run01.nirs');

%% Remove bad channels according SCI
sci = get_SCI(snirf_loaded.data);
disp (sci);
disp(snirf_loaded.data.measurementList.sourceIndex);
disp(snirf_loaded.data.measurementList.detectorIndex);

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
dod_MA_corrected = hmrR_MotionCorrectWavelet(dod, cell(length(dod),1),cell(length(dod),1) ,1.5);

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

%% Reduce extracerebral component with short channels
printf('');