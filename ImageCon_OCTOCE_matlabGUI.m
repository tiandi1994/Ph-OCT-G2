tic

Samplefreq=46816;
Mainfreq=1000; %Vibration Frequency

pixel=1024; %camera pixel
% line=handles.ImageNum; %A-line number
depth=pixel/2; %Image Depth
framenum=512; %frame number
% RefName=sprintf(['','/ref_data']); %reference data name
% coefs = flip([74.6291e-009  -311.7041e-006     1.1153e+000   251.3817e+000]);%from Guan
% coefs = flip([74.6219e-009  -358.7041e-006     1.2253e+000   201.3817e+000]);%from Guan 2
% coefs = flip(pp);
% coefs = flip([79.6291e-009  -355.7041e-006     1.2153e+000   251.3817e+000]);%direct
% coefs = flip([74.6291e-009  -375.7041e-006     1.2553e+000   251.3817e+000]);%Guan refined Untitle3
% coefs = flip([84.5291e-009  -305.7041e-006     1.0553e+000   251.3817e+000]);%Guan unti log only maxim-mean
coefs = flip([84.3319e-009  -304.1041e-006     1.0503e+000 251.3817e+000]); %refine1
% coefs=[-2.475E+1	1.336E+0	-3.767E-4	7.13E-8];%handles
% coefs=handles.coefs;
% coefs=[0 1 0 0];
p=fliplr(coefs); 
K=1:pixel; 
% K=single(K);
KES=polyval(p,-500:1500); %calculate equal step in k space
% KES=single(KES);%Convert double to single data type

KK = interp1(KES,-500:1500,K);
diffMean = (KK(end)-KK(1))/(length(KK)-1);
evenK = diffMean*(cumsum(ones(1,1024))-1)+KK(1);
% RawRef=importdata(RefName)'; %import reference data
% REF=repmat(RawRef,1,line); %generate reference matrix 1024*512

%Parameter for step1
BS_CL=100; %B-scan check line

Deta1= round(Mainfreq*(framenum-2)/Samplefreq)-1; %Vibration Frequency Position
Deta2= round(Mainfreq*(framenum-2)/Samplefreq)+1;
%Parameter for step2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars RawRef
% REF=handles.REF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OCT image construction%%%%%%%%%%%%%%%%%%%%%%%
% load(fileName)
% OCTData = handles.DataCollected;
Spectrum0=permute(single(OCTData),[2 3 1])-repmat(single(REF'),[1 1 framenum]);
CompSpectrum=interp1(KK,Spectrum0,evenK,'linear',0); %Dispersion Compensation
clearvars Spectrum0
temp=fft((CompSpectrum),pixel,1);%fast fourier transform
clearvars filename p CompSpectrum 
% ComplexFrames=temp(depth+1:end,:,:);
ComplexFrames=temp(1:depth,:,:);%dpeth resolved info (Complex matrix 512*512*512)
clearvars temp
figure;imshow(20*log10(abs(squeeze(ComplexFrames(:,:,100)))),[20 80]);title('Structure Image');
colorbar

% ImageBS = mean(abs(ComplexFrames(:,:,:)),3);
ImageBS = abs(ComplexFrames(:,:,1:end-1));

Phframes=angle((ComplexFrames(:,:,2:end-1)).*conj((ComplexFrames(:,:,3:end)))); %Phase image
Phframes = medfilt3(Phframes,[1 1 1]);
% 
% Phframes = angle((ComplexFrames(:,:,2:end-1)))./abs(ComplexFrames(:,:,2:end-1))-...
%     angle((ComplexFrames(:,:,3:end)))./abs(ComplexFrames(:,:,3:end));

% Phframes = -angle((ComplexFrames(:,:,2:end-1)))./ImageBS+...
%     angle((ComplexFrames(:,:,3:end)))./ImageBS;
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NewVibration_ZD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ImCon_Bs=get(handles.Image,'CData'); %Structure image 

A_a1 = Phframes; %Phase Difference
% A_a1=A_a1.*(A_a1<0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b=A_a1<0.5;
% A_a1=A_a1.*b;
% interpl((1:511),A_a1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_1=A_a1>pi;
B_2=A_a1<-pi;
B_3=A_a1~=0;
A_a1=A_a1+2*pi*(B_2-B_1);
A_a1=A_a1.*(B_3-2*B_1).*(B_3-2*B_2);
% A_a1 = medfilt3(A_a1,[3 1 3]);
clearvars B_1 B_2
A_a1 = cumsum(A_a1./1,3);

Ap1=abs(fft(A_a1,[],3));
%&&&&&&&&&&&&&&&& image the vabration &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&        
AmF1=max(Ap1(:,:,Deta1:Deta2),[],3);
MeanAmF=mean(Ap1(:,:,Deta2+150:Deta2+200),3)+0.00001;
Snr=AmF1./MeanAmF;
Snr_TF=Snr>1.25;
AbsImAmpVM1=(AmF1);

% AbsImAmpVM1 = AbsImAmpVM1.*(20*log10(abs(squeeze(ComplexFrames(:,:,100))))>40).*((log10(Snr))>1);


clearvars     Snr_TF

% figure,imshow(ImCon_Bs,[2 4]);%OCT image
% axis on
%%%%%%%%%%%%%%%%%%%%%%%OCE image construction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=3;
strain=zeros(512-m+1,512);

% mask=ImCon_Bs>handles.OCEparameters.StructThresh;
for j=1:512
    for i=1:512-m+1
        strain(i,j)=(6*(2*(i-i)-m+1))*AbsImAmpVM1(i,j)/(m^3-m)*i+(6*(2*(i+1-i)-m+1))...
            *AbsImAmpVM1(i+1,j)/(m^3-m)*(i+1)+(6*(2*(i+2-i)-m+1))*AbsImAmpVM1(i+2,j)/(m^3-m)*(i+2);...
            %+(6*(2*(i+3-i)-m+1))*AbsImAmpVM1(i+3,j)/(m^3-m)*(i+3)+(6*(2*(i+4-i)-m+1))*AbsImAmpVM1(i+4,j)/(m^3-m)*(i+4);
    end
end
Youngs=1./abs(strain);

Y3=193*50*medfilt2(Youngs,[2 2]);
OCEImg=Y3;
mask = (20*log10(abs(squeeze(ComplexFrames(:,:,100))))>40).*(Snr>50);
gra = gradient(AbsImAmpVM1);
% figure;imagesc(medfilt2(abs(1./gra)*80,[3 3]));colormap(jet);title('gradient')
figure;imagesc(Y3,[0 300]);colormap(jet);title('after strain stry');colorbar
figure;imagesc(medfilt2(AbsImAmpVM1,[3 3]),[0 3000]);colormap(jet);title('Vibration Amp')
% a = mean(abs(ComplexFrames(1:510,:,:)),3);
% b = medfilt2(1e4*Youngs./a,[2 2]);
% figure;imagesc(1e4*b,[0 80]);colormap(jet);title('after comp');