%% Notes
%Revised by Duo Zhang, 2018-04-12

%Original file can be found in Labview calibraton Tab.

%This file generates coefs for dispersion calibration.

%% Parameters setting
clc
pixel=1024;
%rem_stt and rem_end set the region for credible peaks
rem_stt = 250; %Start points for credible peaks
rem_end =50;    %End points for credible peaks
height = 2;
%Open data file
fileID = fopen('2.txt','r');
A = fscanf(fileID,'%u');
B = reshape(A,[2,1024]);
%Data stored in sp
sp = B(2,:);
fftsp0 = log10(abs(fft(sp)));%Original fft
sp([1:rem_stt,pixel-rem_end+1:end])=0;
[~,locs] = findpeaks(sp,'MinPeakProminence',35);
figure;findpeaks(sp,'MinPeakProminence',35);
%% Coefs finding
%Note: Peak numbers are not integer values. For polynomial fitting,
%continous values rather than discrete values are recommended. 
newx = linspace(locs(1),locs(end),locs(end)-locs(1)+1);%location of the peaks
figure;plot(locs,1:length(locs),'o');%Location of the peaks & peak No. before interp

pkNo = interp1(locs,1:length(locs),newx,'spline');%corresponded peak number
hold on;plot(newx,pkNo,':.') %Location of the peaks & peak No. after interp
xlabel('peak location')
ylabel('Peak No.')
%Convert space
relatedWaveL = (locs(end)-locs(1))/(pkNo(end)-pkNo(1))...
    *(pkNo-pkNo(1))+pkNo(1);%Location related peak No.
pp = polyfit(relatedWaveL,newx,3);%peak number & peak location fitting
ngg = polyval(pp,1:1024);
SS=interp1(1:1024,B(2,:),ngg,'linear',0);
II=log10(abs((fft(SS))));
figure;plot(II);hold on;plot(fftsp0);legend('after Comp','before Comp')
ylim([1 6])
stringDisp = sprintf(['\n','coefs = ',num2str(pp),'\n\n','Remember to filp it for OCT/OCE processing when in use','\n']);
fprintf(stringDisp)
