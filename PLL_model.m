% mmWave PLL models

DF = 100e3; % SCS
CBW = 115e6; % BW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   29.55 GHz PLL model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSD0 = 1585; % 32 dB
f_zn = [3e3 550e3 280e6];
a_zn = [2.37 2.7 2.53];
f_pm = [1 1.6e6 30e6];
a_pm = [3.3 3.3 1];

f0 = 100:DF:CBW;
Lenf0 = length(f0);
S_f0 = zeros(size(Lenf));

for ii = 1:Lenf0
    S_f0
