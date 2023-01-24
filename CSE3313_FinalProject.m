%Name: Athiya Manoj
%UTA ID: 1001843726

clc
clear all
close all

%asking for user input
userInput = "Please enter the file (.wav)";
filename = input(userInput,"s");

%read the audio file and get sampling frequency 
[y,Fs] = audioread(filename);

%%Reading the file
%Length of signal
L = length(y);
%Sampling period
T = 1/Fs;

%Filtering the audio in time domain
figure(1),spectrogram(y, 1024,[], 1024, Fs), title('Spectrogram of Input Signal');

%Plotting DFT to find frequency
f = linspace(-Fs/2, Fs/2, L);
figure(2)
plot(f,fftshift(abs(fft(y))))
title('DFT of Input audio')
xlabel('Frequency in Hz')
grid on

%conversation is at about 0.016 with harmonic at 318 Hz
%the high freq tone is at 4353 Hz
wp = 300;   %will keep conversation and harmonic
ws = 3500;  %will eliminate beep at 4kHz
Rp = 3;
Rs = 50;


%gives normalized frequencies
Wp = wp/(Fs/2);
Ws = ws/(Fs/2);

[N,Wc] = buttord(Wp,Ws,Rp,Rs);

%Filter needs an order of 2 and Wc = 0.0557
%Checking actual frequency
wc = Wc*(Fs/2);

%Filter Design
[b, a] = butter(N, Wc, 'low');

%Viewing the filter response from 0 to pi for 1000 points
[H,W] = freqz(b,a,1000);

%digital frequency linear plots
figure(3), plot(W, abs(H)), title('Digital Filter frequency response'), xlabel('Frequency in rad/sec');
figure (4), plot(W, 20*log10(abs(H))), title('Digital Filter frequency response'), xlabel('Frequency in rad/sec');

%continuous freq
f1 = linspace(0,Fs/2,1000);
figure(5),plot(f1,20*log10(abs(H))), title('Actual Continuous frequency'),xlabel('Frequency in rad/sec'), ylabel('Gain in dB');
 
output = filter(b,a,y);

figure(6), spectrogram(output, 1024, [], 1024, Fs), title(' Spectrogram of output using filter');


%implement in a computer with difference equations
na = length(a); %denominator coefficients
nb = length(b); %numerator coefficient
y1 = zeros(1,L);
for i = 1:L %step over the entire input file 
  for j = 2:na %step over denominator 'a' coefficients
      yindex = i-j;
      if yindex <= 0 %ignore attempts to get -ve indexes from y
          y1(i) = y1(i)+0;
      else
          y1(i) = y1(i) - a(j)*y1(yindex+1);
      end
  end 
  for j = 1:nb %step over denominator 'b' coefficients
      xindex = i-j;
      if xindex <= 0 %ignore attempts to get -ve indexes from x
          y1(i) = y1(i)+0;
      else
          y1(i) = y1(i) + b(j)*y(xindex+1);
      end
  end
end

figure(7), plot(abs(fftshift(fft(y1)))), title('DFT of output'), xlabel('Frequency in Hz'), grid on;

figure(8),spectrogram(y1,1024,[],1024,Fs), title('Spectrogram of output using Difference Equations');

figure(9), subplot(3,1,1),plot(y),title('Input');
figure(9), subplot(3,1,2),plot(output),title('Output using filter() command');
figure(9), subplot(3,1,3),plot(y1),title('Output using Difference Equation');

%listening to filtered audio
sound(output,Fs)
audiowrite('filteredconversation.wav',output, Fs)
