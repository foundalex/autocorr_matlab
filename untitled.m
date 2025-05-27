clear all;
sintablen = 4096;
for j = 0:sintablen-1
    COSTAB_FOR(j+1) = cos(2*pi*j/sintablen);
end

SINTAB = (-sin(2*pi*(0:sintablen-1)./sintablen));
COSTAB = (cos(2*pi*(0:sintablen-1)./sintablen));

% n = 1:10000;
% time_base=0:10000-1;
% correction_signal=exp(i*(0.0006283)*time_base);
% correction_signal1=exp(i*(0.00006283)*time_base);
% plot(n, real(correction_signal), n, real(correction_signal1))

 error_costab = COSTAB - COSTAB_FOR;

fs1 = 20000000;
F_required = 100;
index = 1; 
step = ((F_required/fs1)*sintablen);
phaser = 0;
for i = 1:40000
    dds_cos1(i) = COSTAB(round(index));
    dds_sin1(i) = SINTAB(round(index));

    index = index+step;

    if index>sintablen
      index = index - sintablen;
    end
    indexii(i) = (index);

%     sig(i) = cordicrotate(phaser, 1, 14);
%     phaser = phaser + 0.6283;
%         if (phaser < -pi)
%             phaser = phaser + 2*pi;
%         elseif (phaser > pi)
%             phaser = phaser - 2*pi;
%         end
end
% 
figure(11)
n = (1:40000).';
% 
%     spectrumScope = spectrumAnalyzer(SampleRate=1000, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([sig.']);

plot(n, dds_cos1)
% plot(n, indexii*2*pi/1/fs1)
% title('Адрес таблицы от фазы')
% xlabel('Номер отсчета') 
% ylabel('Адрес из таблицы') 