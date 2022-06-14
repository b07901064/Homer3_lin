function output = get_SCI(data, nb_channels)
%get_SCI Summary of this function goes here
%   Detailed explanation goes here
%   filtered both photodetected signals between 0.5 and 2.5 Hz \
%   to preserve only the cardiac component
%   and normalized the resulting signals to balance any difference between their amplitude.
%   computed the cross-correlation and we extracted the value at a time lag of 0

sci = zeros(nb_channels,1);
data_bpf = hmrR_BandpassFilt(data, 0.5 , 2.5);
data_bpf = normalize(data_bpf.dataTimeSeries);
for i = 1:nb_channels
    [r,lags] = xcorr(data_bpf(:,i),data_bpf(:,i+ nb_channels), 'coeff');
    sci(i) = r(lags == 0);
end
output = sci;
end

