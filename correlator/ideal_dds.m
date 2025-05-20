
% function [s]=ideal_dds(fs,bp,bd,fd,ps,N)
% N : length of time domain signal (points)
% bp: number of bits (precision) of sine phase lookup table (bits)
% bd: DAC resolution (bits)
% fs: sampling frequency (Hz)
% fd: desired output frequency (Hz)
% ps: desired starting phase (deg)
%
%
N = 3221;
fd = 100000;
fs = 20000000;
bp = 12;
bd = 16;
ps = 90;
% calculate frequency tuning word (FTW) from desired frequency
ftw=round(fd/fs*2^bp);

% calculate actual synthesized frequency (Hz)
f0=ftw*fs/2^bp;

% calculate phase tuning word (PTW) for start phase
ptw=round(mod(ps/360*2^bp,2^bp));

% calculate actual synthesized phase start (scaled rad)
p0=ptw/2^bp;

% create time vector for sinusoid phase (s)
t=1/fs*[0:N-1];

% phase of sinusoid (rad)
phase=round(mod(f0*t+p0,1)*2^bp)/2^bp*2*pi;

% quantized sinusoid value
s=round(sin(phase)*(2^bd-1));

% normailzed sinusoid
s=s/(2^bd-1);

%     spectrumScope = spectrumAnalyzer(SampleRate=fs, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([s.']);
% fprintf(‘Synthesized frequency = %0.4f kHz at starting angle %0.4f degrees\n’,f0/1e3,p0*360);

% end
