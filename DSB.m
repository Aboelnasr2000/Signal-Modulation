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
signalSCFD=real(fftshift(fft(signalSC)));
n=length(signalSCFD);
xNew=linspace(-newFS/2, newFS/2,length(signalSCFD));
figure(2);
subplot(2,1,1);
plot(xNew,abs(signalSCFD),'r');
ylim([0 4000]);
xlabel('Frequency');
title('DSB-SC Frequency domain');

%Max amplitude of filtered audio
Am=max(abs(YtimeDomain));
Ac=2*Am;

%DSB-TC modulation
signalTC=(Ac + YtimeDomain).*cos(2*pi*fc*t); 

%Sketch the modulated DSB-TC in frequency domain
signalTCFD=real(fftshift(fft(signalTC)));
subplot(2,1,2);
plot(xNew,abs(signalTCFD),'r');
ylim([0 5000]);
xlabel('Frequency');
title('DSB-TC Frequency domain');

%Envelope detector for both modulation types 
%Envelope can't be used with DSB-SC
envelopeSC=abs(hilbert(real(signalSC)));
envelopeTC=abs(hilbert(real(signalTC)));

%DC bias
% SCdemodulated = detrend(envelopeSC); 
TCdemodulated = detrend(envelopeTC); 

%Resample
SCdemodulated =resample(envelopeSC,Fs,5*fc);
TCdemodulated=resample(TCdemodulated,Fs,5*fc);

%Play demodulated signal
%Save Sound Data
play_DemodEnvY=real(TCdemodulated);
play_DemodEnvFs=Fs;

%For matrix dimensions 
time = transpose(time);
time = [0; time];

figure(3);
subplot(2,1,1);
plot(time,SCdemodulated);
xlabel('Time');
title('Envelope Detection of DSB-SC (Invalid)');

subplot(2,1,2);
plot(time,TCdemodulated);
xlabel('Time');
title('Envelope Detection of DSB-TC ');

%Coherent detection DSB-SC
%SNR=0
%Add noise 0 dB
SC_SNR0 = awgn(signalSC,0);
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
figure(5);
subplot(2,1,1);
plot(time,real(coherentFilteredSNR0));   
xlabel('Time');
title('Coherent detector DSB-SC with SNR=0');
x = [0 x];
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSNR0)))),'r');                     
xlabel('Frequency');
title('Coherent detector Spectrum DSB-SC with SNR=0');

% SNR=10 
SC_SNR10 = awgn(signalSC,10);
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
figure(6);
subplot(2,1,1);
plot(time,real(coherentFilteredSNR10));
xlabel('Time');
title('Coherent detector DSB-SC with SNR=10');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSNR10)))),'r');
xlabel('Frequency');
title('Coherent detector Spectrum DSB-SC with SNR=10')

%SNR=30
SC_SNR30 = awgn(signalSC,30);
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
figure(7);
subplot(2,1,1);
plot(time,real(coherentFilteredSNR30));
xlabel('Time');
title('Coherent detector DSB-SC with SNR=30');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentFilteredSNR30)))),'r');
xlabel('Frequency');
title('Coherent detector Spectrum DSB-SC with SNR=30');

%Carrier Frequency Offset
%Coherent detection with frequency error F=100.1 KHz
fcerror=fc*1.001;
coherentError=signalSC.*cos(2*pi*fcerror*t);
%Remove all frequencies greater than 4000 Hz
%Zero-phase digital filtering
coherentErrorFiltered = filtfilt (b, a, coherentError).*2;

%resample to original signal
coherentErrorFiltered =resample(coherentErrorFiltered,Num,Den);

%Save Sound Data
play_DemodCohFreqErrorY=real(coherentErrorFiltered);
play_DemodCohFreqErrorFs=Fs;

%Plot Time & Frequency Domain
figure(8);
subplot(2,1,1);
plot(time,real(coherentErrorFiltered));
xlabel('Time');
title('Coherent detection with frequency error');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentErrorFiltered)))),'r');
xlabel('Frequency');
title('Coherent detection Spectrum with frequency error')

%Coherent detection with phase error = 20
%Convert to radian
PhaseError=(20*2*pi)/180;
coherentPhaseError=signalSC.*cos(2*pi*fc*t+PhaseError);
%Remove all frequencies greater than 4000 Hz
%Zero-phase digital filtering
coherentPhaseErrorFiltered = filtfilt (b, a, coherentPhaseError).*2; 

%Decrease the sampling frequency again, (Return to original fs)
coherentPhaseErrorFiltered =resample(coherentPhaseErrorFiltered,Num,Den);

%Save Sound Data
play_DemodCohPhaseErrorY=real(coherentPhaseErrorFiltered);
play_DemodCohPhaseErrorFs=Fs;


%Plot Time & Frequency Domain
figure(9);
subplot(2,1,1);
plot(time,real(coherentPhaseErrorFiltered));    
xlabel('Time');
title('Coherent detection with phase error 20');
subplot(2,1,2);
plot(x,abs(real(fftshift(fft(coherentPhaseErrorFiltered)))),'r');    
xlabel('Frequency');
title('Coherent detection Spectrum with phase error 20')

%Menu to play sound
flag=1;
while flag==1
    choice = menu('Play Audio','Original','Filtered','Envelope DSB-SC','Envelope DSB-TC','Coherent DSB-SC SNR=0','Coherent DSB-SC SNR=10','Coherent DSB-SC SNR=30','Coherent DSB-SC Frequency Error=0.1kHz','Coherent DSB-SC Phase Error=20','Close');
    if choice==1
        sound(play_originalY,play_originalFs);
    else if choice==2
        sound(play_FilteredY,play_FilteredFs); 
    else if choice==3
        msgbox('Envelope Detector invalid with DSB-SC', 'Invalid','error');
    else if choice==4
        sound(play_DemodEnvY,play_DemodEnvFs);
    else if choice==5
        sound(play_DemodCohSNR0Y,play_DemodCohSNR0Fs);
    else if choice==6
        sound(play_DemodCohSNR10Y,play_DemodCohSNR10Fs);
    else if choice==7
        sound(play_DemodCohSNR30Y,play_DemodCohSNR30Fs);
    else if choice==8
        sound(play_DemodCohFreqErrorY,play_DemodCohFreqErrorFs);
    else if choice==9
        sound(play_DemodCohPhaseErrorY,play_DemodCohPhaseErrorFs);
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
end
        
    