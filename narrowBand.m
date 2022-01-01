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
[Num,Den]=rat(newFS/Fs);
YtimeDomain=resample(YtimeDomain,Num,Den); 

t= linspace(0,length(YtimeDomain)/newFS,length(YtimeDomain));
x = linspace(-newFS/2,newFS/2,length(YtimeDomain));
t=transpose(t); 

%integration of m(t)
m_integration=cumsum(YtimeDomain)/newFS;
%delta fmax
Am=max(abs(YtimeDomain));
%Assume beta is 0.02 <<1 ---> NBFM
%Condition for NBFM beta<<1
beta = 0.02;
%Calculate Kf
kf=(beta*4000)/Am;
NBFM=cos(2*pi*fc*t + kf*2*pi*m_integration);

%Sketch the modulated NBFM frequency domain
figure(2)
subplot(2,1,1);
plot(x,abs(real(fftshift(fft(NBFM)))),'r');
ylim([0 3000]);
xlabel('Frequency');
title('NBFM Frequency Domain');
output=diff(NBFM)*250;
envelope=abs(hilbert(real(output)));

%DC bias
demodulated=detrend(envelope);

%Resample to original
[Num,Den]=rat(Fs/newFS);
demodulated=resample(demodulated,Num,Den);
sound(real(demodulated),Fs);

subplot(2,1,2);
plot(time,demodulated);
ylim([-0.3 0.3]);
xlabel('time');
title('Demodulated NBFM Time Domain');

