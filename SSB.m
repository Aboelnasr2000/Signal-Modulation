clear;
clc;
close all;
%stop any playing audio
clear sound;

%Get Audio File 
audio = uigetfile;   
%y --> audio, fs --> sampling frequency
[y,Fs] = audioread(char(audio));

%Save Sound Data
play_originalY=real(y);
play_originalFs=Fs;

%get fourier transform of audio signal
Y=fft(y); 
%Center signal at zero
Yshift=fftshift(Y);

%Plot fourier transformed signal
x = linspace( -Fs/2, Fs/2,length(Yshift));
figure(1);
subplot(3,1,1);
plot(x,abs(Yshift),'r');
xlabel('Frequency');
title('Original Spectrum');

%Ideal Filter
Yshift(x<-4000)=0;
Yshift(x>4000)=0;

subplot(3,1,2);
plot(x,abs(real(Yshift)),'r');
xlabel('Frequency');
title('Filtered Spectrum');

%frequency domain to time domain
%unshift the signal
Yunshift=ifftshift(Yshift);
%inverse fourier transfrom
YtimeDomain=ifft(Yunshift); 

time = linspace(0,length(YtimeDomain)/Fs, length(YtimeDomain));
subplot(3,1,3);
plot(time,real(YtimeDomain));
xlabel('Time');
title('Filtered Audio Time Domain');

%Play Filtered Audio(real part)
%Save Sound Data
play_FilteredY=real(YtimeDomain);
play_FilteredFs=Fs;

fc=100000;

newFS=5*fc;
%Get simplest rational fraction
[Numerator,Denominator]=rat(newFS/Fs);
YtimeDomain=resample(YtimeDomain,Numerator,Denominator); 

t= linspace( 0,  length(YtimeDomain)/newFS,  length(YtimeDomain) );
%For matrix dimensions 
t=transpose(t);

%DSB-SC modulation
signalSC=YtimeDomain.*cos(2*pi*fc*t); 

%Sketch the modulated DSB-SC in frequency domain
signalSCFD=fftshift(fft(signalSC));
xNew=linspace(-newFS/2, newFS/2,length(signalSCFD));
figure(2);
subplot(2,1,1);
plot(xNew,abs(signalSCFD),'r');
ylim([0 4000]);
xlabel('Frequency');
title('DSB-SC Frequency domain');

%ideal filter
signalSCFD(xNew>100000)=0;
signalSCFD(xNew<-100000)=0;

%Sketch the modulated SSB-SC in frequency domain
subplot(2,1,2);
plot(xNew,abs(real(signalSCFD)),'r');
ylim([0 4000]);
xlabel('Frequency');
title('SSB-SC Frequency domain');

time = transpose(time);
time = [0; time];

%frequency domain to time domain
%unshift the signal
signalSCFD_Unshift=ifftshift(signalSCFD);
%inverse fourier transfrom
signalSCTD=ifft(signalSCFD_Unshift); 

%Coherent detection SSB-SC
coherent_SSB_Ideal=signalSCTD.*cos(2*pi*fc*t);

%Get b & a of filter using butterworth order=5 & multiply by 2 due to
%division by 2 from fourier transfrom of cos
[b, a] = butter (5, 4000.*2./newFS);
%Zero-phase digital filter & multiply by 2
coherentFilteredSSB_Ideal = filtfilt (b, a, coherent_SSB_Ideal).*2; 

%Resample to original frequency
[Num,Den]=rat(Fs/newFS);
coherentFilteredSSB_Ideal =resample(coherentFilteredSSB_Ideal,Num,Den);

%Save Sound Data
play_DemodIdealSCY=real(coherentFilteredSSB_Ideal);
play_DemodIdealSCFs=Fs;

%Plot Time & Frequency Domain
figure(3);
subplot(2,1,1);
plot(time,real(coherentFilteredSSB_Ideal));   
xlabel('Time');
title('Coherent detector SSB-SC Ideal Filter with no noise');
x = [0 x];
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSSB_Ideal)))),'r');                     
xlabel('Frequency');
title('Coherent detector Spectrum SSB-SC Ideal Filter with no noise');


%butterworth 4th order filter
w=[(100000-4000)*2/newFS 100000*2/newFS];
[b, a] = butter (4, w);
signal_butter = fftshift(fft(signalSC));
signalSC_Butter = filtfilt (b, a, signalSC).*2; 
% signalSC_Butter=filter(b,a,signalSC);

%Sketch the modulated SSB-SC in frequency domain
figure(4);
subplot(2,1,1);
plot(xNew,abs(signal_butter),'r');
ylim([0 4000]);
xlabel('Frequency');
title('DSB-SC Frequency domain');
subplot(2,1,2);
xNew=linspace(-newFS/2, newFS/2,length(signalSC_Butter));
plot(xNew,abs(real(fftshift(fft(signalSC_Butter)))),'r');
ylim([0 4000]);
xlabel('Frequency');
title('SSB-SC Butterworth Filtered Frequency domain');

%frequency domain to time domain
%unshift the signal
signalSCFD_Unshift=ifftshift(signalSCFD);
%inverse fourier transfrom
signalSCTD=ifft(signalSCFD_Unshift); 

%Coherent detection SSB-SC
coherent_SSB_Butter=signalSCTD.*cos(2*pi*fc*t);

%Get b & a of filter using butterworth order=5 & multiply by 2 due to
%division by 2 from fourier transfrom of cos
[b, a] = butter (5, 4000.*2./newFS);
%Zero-phase digital filter & multiply by 2
coherentFilteredSSB_Butter = filtfilt (b, a, coherent_SSB_Butter).*2; 

%Resample to original frequency
[Num,Den]=rat(Fs/newFS);
coherentFilteredSSB_Butter =resample(coherentFilteredSSB_Butter,Num,Den);

%Save Sound Data
play_DemodButterY=real(coherentFilteredSSB_Butter);
play_DemodButterFs=Fs;

%Plot Time & Frequency Domain
figure(5);
subplot(2,1,1);
plot(time,real(coherentFilteredSSB_Butter));   
xlabel('Time');
title('Coherent detector SSB-SC Butterworth Filtered');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSSB_Butter)))),'r');                     
xlabel('Frequency');
title('Coherent detector Spectrum SSB-SC Butterworth Filtered');

%Coherent Detection SSB-SC
%SNR=0
%Add noise 0 dB
signalSSBSCTD=ifft(ifftshift(signalSCFD));
SC_SNR0 = awgn(signalSSBSCTD,0);
coherent_SNR0=SC_SNR0.*cos(2*pi*fc*t);

%Get b & a of filter using butterworth order=5 & multiply by 2 due to
%division by 2 from fourier transfrom of cos
[b, a] = butter (5, 4000.*2./newFS);
%Zero-phase digital filter & multiply by 2
coherentFilteredSNR0 = filtfilt (b, a, coherent_SNR0).*2; 

%Resample to original frequency
[Num,Den]=rat(Fs/newFS);
coherentFilteredSNR0 =resample(coherentFilteredSNR0,Num,Den);

%Save Sound Data
play_DemodCohSNR0Y=real(coherentFilteredSNR0);
play_DemodCohSNR0Fs=Fs;

%Plot Time & Frequency Domain
figure(6);
subplot(2,1,1);
plot(time,real(coherentFilteredSNR0));   
xlabel('Time');
title('Coherent detector SSB-SC with SNR=0');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSNR0)))),'r');                     
xlabel('Frequency');
title('Coherent detector Spectrum SSB-SC with SNR=0');

% SNR=10 
SC_SNR10 = awgn(signalSSBSCTD,10);
coherent_SNR10=SC_SNR10.*cos(2*pi*fc*t);
%Remove all frequencies greater than 4000 Hz
%Zero-phase digital filtering
coherentFilteredSNR10 = filtfilt (b, a, coherent_SNR10).*2;                   

%Decrease the sampling frequency again, (Return to original fs)
[Num,Den]=rat(Fs/newFS);
coherentFilteredSNR10 =resample(coherentFilteredSNR10,Num,Den);

%Save Sound Data
play_DemodCohSNR10Y=real(coherentFilteredSNR10);
play_DemodCohSNR10Fs=Fs;

%Plot Time & Frequency Domain
figure(7);
subplot(2,1,1);
plot(time,real(coherentFilteredSNR10));
xlabel('Time');
title('Coherent detector SSB-SC with SNR=10');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSNR10)))),'r');
xlabel('Frequency');
title('Coherent detector Spectrum SSB-SC with SNR=10')

%SNR=30
SC_SNR30 = awgn(signalSSBSCTD,30);
coherent_SNR30=SC_SNR30.*cos(2*pi*fc*t);
%Remove all frequencies greater than 4000 Hz
%Zero-phase digital filtering
coherentFilteredSNR30 = filtfilt (b, a, coherent_SNR30).*2; 
%Decrease the sampling frequency again, (Return to original fs)
coherentFilteredSNR30 =resample(coherentFilteredSNR30,Num,Den);

%Save Sound Data
play_DemodCohSNR30Y=real(coherentFilteredSNR30);
play_DemodCohSNR30Fs=Fs;

%Plot Time & Frequency Domain
figure(8);
subplot(2,1,1);
plot(time,real(coherentFilteredSNR30));
xlabel('Time');
title('Coherent detector SSB-SC with SNR=30');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSNR30)))),'r');
xlabel('Frequency');
title('Coherent detector Spectrum SSB-SC with SNR=30');


%%SSB-TC
%Max amplitude of filtered audio
Am=max(abs(YtimeDomain));
Ac=2*Am;

%DSB-TC modulation
signalTC=(Ac + YtimeDomain).*cos(2*pi*fc*t); 

%Sketch the modulated DSB-TC in frequency domain
signalTCFD=fftshift(fft(signalTC));
xNew=linspace(-newFS/2, newFS/2,length(signalTCFD));
figure(9);
subplot(2,1,1);
plot(xNew,abs(signalTCFD),'r');
ylim([0 5000]);
xlabel('Frequency');
title('DSB-TC Frequency domain');

%ideal filter
signalTCFD(xNew>100001)=0;
signalTCFD(xNew<-100001)=0;

%Sketch the modulated SSB-SC in frequency domain
subplot(2,1,2);
xNew=linspace(-newFS/2, newFS/2,length(signalTCFD));
plot(xNew,abs(real(signalTCFD)),'r');
ylim([0 4000]);
xlabel('Frequency');
title('SSB-TC Frequency domain');

%Get time domain
signalTCTD=ifft(ifftshift(signalTCFD));

%Envelope detector
envelopeTC=abs(hilbert(real(signalTCTD)));

%DC Bias
TCdemodulated = detrend(envelopeTC);

%Resample
TCdemodulated=resample(TCdemodulated,Fs,5*fc);

%Play demodulated signal
%Save Sound Data
play_DemodCohTCY=real(TCdemodulated);
play_DemodCohTCFs=Fs;

figure(10);
subplot(2,1,1);
plot(time,TCdemodulated);
xlabel('Time');
title('Envelope Detection of SSB-TC');

%Menu to play sound
flag=1;
while flag==1
    choice = menu('Play Audio','Original','Filtered','Coherent Detection Ideal Filter SSB-SC','Coherent Detection Butterworth Filter SSB-SC','Coherent SSB-SC SNR=0','Coherent SSB-SC SNR=10','Coherent SSB-SC SNR=30','Coherent Ideal Filter SSB-TC','Close');
    if choice==1
        sound(play_originalY,play_originalFs);
    else if choice==2
        sound(play_FilteredY,play_FilteredFs); 
    else if choice==3
        sound(play_DemodIdealSCY,play_DemodIdealSCFs);
    else if choice==4
        sound(play_DemodButterY,play_DemodButterFs);
    else if choice==5
        sound(play_DemodCohSNR0Y,play_DemodCohSNR0Fs);
    else if choice==6
        sound(play_DemodCohSNR10Y,play_DemodCohSNR10Fs);
    else if choice==7
        sound(play_DemodCohSNR30Y,play_DemodCohSNR30Fs);
    else if choice==8
        sound(play_DemodCohTCY,play_DemodCohTCFs);
        else
            flag=0;
        end
        end
        end
        end
        end
        end
        end
    end
end
        
    






















