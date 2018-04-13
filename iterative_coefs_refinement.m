%% Notation
%File revised by Duo Zhang, 2018-04-13. 

%Original file is programmed by Guangying Guan, with the name as
%'SearchOptimizeParaNoguiyimm_FWHM_1'

%This program can refine coefs generated by coefsGeneration.m

%Set loop valuse W1 W2 W3 respectivly which correspond to coef(1:3) in a
%resonable range.

%% parameters setting
pixel=1024;

W0=251.3817e+000; %Coefs(4) which is not used for refinement 
FWHMFactor=1.0556e+000; % Factor used to find width of peak (not in use in this version)
rrr=3;%Number for determining continuity of conotonicity for peak

%Open file for coefs refinement
fileID = fopen('2.txt','r');
A = fscanf(fileID,'%u');
B = reshape(A,[2,1024]);
specO=B(2,:);

P1=550; %P1 and P2 set borders to find location of peak
P2=700;
M=300; %Set initial value for FWHM which is larger than normal (but this property is not in use in this version)
OPMW1=0; %Zero coefs at the beginning
OPMW2=0;
OPMW3=0;
MM=0; %Set intial value for maimum peak height
%% Iterative coefs finding
for W1=74.6291e-009:0.1e-9:84.6291e-9
    W1;
    for W2=-405.7041e-006:1e-6:-305.7041e-006
        W22=W2;   
        
        for W3=0.7153e+000:0.01:1.7153
            W33=W3;
            for ll=1:1
                
                
                pp2 = [W1 W22 W33 W0];
                
                waveL = polyval(pp2,-500:1500); %waveLength
                
                KK = interp1(waveL,-500:1500,1:pixel);
                diffMean = (KK(end)-KK(1))/(length(KK)-1);
                evenK = diffMean*(cumsum(ones(1,pixel))-1)+KK(1);
                
try
                SpecCI=interp1(KK,specO,evenK,'spline');
catch
    break
end

                IC=20*log10(abs(fftshift(fft(SpecCI))));

                [M1,IM1]=max(IC(P1:P2));

                 ss=0;
                for jj=P1+IM1-rrr:P1+IM1-1
                    if IC(jj)-IC(jj-1)<0
                       ss=1;
                    end
                end
                
                if ss==1 
                    break;
                end
                
                for jj=P1+IM1+1:P1+IM1+rrr
                    if IC(jj)-IC(jj-1)>0
                       ss=1;
                    end
                end
                
                if ss==1 
                    break;
                end

               

                iii=IM1;
                while(IC(P1+iii-1)>M1/FWHMFactor)
                    iii=iii+1;
                end
                IM1B=iii;
                iii=IM1;
                while(IC(P1+iii-1)>M1/FWHMFactor)
                    iii=iii-1;
                end
                IM1S=iii;

                M1=M1-mean(IC(630:700));
             
                SUM_FWHM=IM1B-IM1S;
                if (M1>MM)
                    MM=M1;
                    OPMW1=W1;
                    OPMW2=W22;
                    OPMW3=W33;
                    M=SUM_FWHM;
                    figure(30)
                    plot((IC(:)));
                    title(num2str(pp2))
                    drawnow
                    
                end
            end
        end
    end
end
pp2=[OPMW1 OPMW2 OPMW3 W0]%result




