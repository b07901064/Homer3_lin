function output = get_SCI(data)
%REMOVECH_SCI Summary of this function goes here
%   Detailed explanation goes here
sci = zeros(9,1);
for i = 1:9
    [r,lags] = xcorr(data.dataTimeSeries(:,i),data.dataTimeSeries(:,i+9));
    sci(i) = r(lags == 0);
end
output = sci;
end

