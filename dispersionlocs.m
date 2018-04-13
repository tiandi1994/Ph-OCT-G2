%%
pixel=1024;
rem_stt = 250;
rem_end =50;
height = 2;

fileID = fopen('2.txt','r');
A = fscanf(fileID,'%u');
B = reshape(A,[2,1024]);

sp = B(2,:);
fftsp0 = log10(abs(fft(sp)));
sp([1:rem_stt,pixel-rem_end+1:end])=0;
[~,locs] = findpeaks(sp,'MinPeakProminence',35);
figure;findpeaks(sp,'MinPeakProminence',35);
%%
%Note: Peak numbers are not integer values. For polynomial fitting,
%continous values rather than discrete values are recommended. 
newx = linspace(locs(1),locs(end),locs(end)-locs(1)+1);%location of the peaks
figure;plot(locs,1:length(locs),'o');%Location of the peaks & peak No. before interp

pkNo = interp1(locs,1:length(locs),newx,'spline');%corresponded peak number
hold on;plot(newx,pkNo,':.') %Location of the peaks & peak No. after interp
xlabel('peak location')
ylabel('Peak No.')
%Convert 
relatedWaveL = (locs(end)-locs(1))/(pkNo(end)-pkNo(1))...
    *(pkNo-pkNo(1))+pkNo(1);%Location related peak No.
pp = polyfit(relatedWaveL,newx,3);%peak number & peak location fitting
% pp = [97.6097e-009  -349.5812e-006     1.1438e+000   205.5863e+000];%result1
% pp = [92.6097e-009  -344.5812e-006     1.1538e+000   205.5863e+000];%result2
ngg = polyval(pp,1:1024);
SS=interp1(1:1024,B(2,:),ngg,'linear',0);
II=log10(abs((fft(SS))));
figure;plot(II);hold on;plot(fftsp0);legend('after Comp','before Comp')
ylim([1 6])
