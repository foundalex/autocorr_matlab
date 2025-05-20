function [dds_cos, dds_sin] = dds_int(angle_rad_int, angle_in_rad_int_16, n)
% 1/4 period of sine
lines = readlines("test_signals\rot_lut.txt");
% a = importdata("test_signals\rotate_mult_i.txt");

% convert int sine to double
for i = 1:512
    str = num2str(lines(i));
    r = str - '0';
    rr = dec2bin(r);
    str_q = [rr(17), rr(18), rr(19), rr(20), rr(21), rr(22), rr(23), rr(24), rr(25), rr(26), rr(27), rr(28), rr(29), rr(30), rr(31), rr(32)];
    str_i = [rr(1), rr(2), rr(3), rr(4), rr(5), rr(6), rr(7), rr(8), rr(9), rr(10), rr(11), rr(12), rr(13), rr(14), rr(15), rr(16)];

    str_i_array(i) = bin2dec(str_i);
    str_q_array(i) = bin2dec(str_q);
end
%%
% cfgHE = wlanHESUConfig;
% cbw = cfgHE.ChannelBandwidth;
% waveform = wlanWaveformGenerator(1,cfgHE);
% %Получите индексы поля WLAN, которые содержат модулируемые поля L-SIG и RL-SIG.
% 
% ind = wlanFieldIndices(cfgHE);
% rxLSIG = waveform(ind.LSIG(1):ind.RLSIG(2),:);
% %Выполните демодуляцию ортогонального мультиплексирования деления частоты (OFDM), чтобы извлечь поле L-SIG.
% 
% lsigDemod = wlanHEDemodulate(rxLSIG,'L-SIG',cbw);
% %Насчитайте символы L-SIG и RL-SIG, возвратите предHE информация о OFDM и извлеките демодулируемые символы L-SIG.
% 
% lsigDemodAverage = mean(lsigDemod,2);
% preHEInfo = wlanHEOFDMInfo('L-SIG',cbw);
% lsig = lsigDemodAverage(preHEInfo.DataIndices,:);
% %Восстановите биты информации о L-SIG и другую информацию, не приняв шума канала. Отобразите результат проверки четности.
% 
% noiseVarEst = 0;
% [bits,failCheck,info] = wlanLSIGBitRecover(lsig,noiseVarEst);
% disp(failCheck);
%%
pi_int = 1608;
% pi_int = 180*512;
%
index = 0;
next_phase_correction = 0;

for i = 1:n

    if (angle_rad_int < 0)
        next_phase_correction = next_phase_correction - round((angle_in_rad_int_16*512),3);
%             next_phase_correction = next_phase_correction - (angle_in_rad_int_16/512);
        if (next_phase_correction < -pi_int)
            next_phase_correction = next_phase_correction + 2*pi_int;
        elseif (next_phase_correction > pi_int)
            next_phase_correction = next_phase_correction - 2*pi_int;
        end
    end

    if abs(next_phase_correction) <= pi_int/4
        if (next_phase_correction <= 0)
            quadrant = 4;
        else
            quadrant = 0;
        end
        index = floor(abs(next_phase_correction));
    elseif abs(next_phase_correction) <= pi_int/2
        if (next_phase_correction < 0)
            quadrant = 5;
        else
            quadrant = 1;
        end
        index = floor(pi_int/2 - abs(next_phase_correction));
    elseif abs(next_phase_correction) <= pi_int*3/4
        if (next_phase_correction < 0)
            quadrant = 6;
        else
            quadrant = 2;
        end
          index = floor(abs(next_phase_correction) - pi_int/2);
    else
        if (next_phase_correction < 0)
            quadrant = 7;
        else
            quadrant = 3;
        end
        index = round(pi_int - abs(next_phase_correction));
    end
    
    sin1(i) = str_i_array(index+1);
    cos1(i) = str_q_array(index+1);

        if (quadrant == 0)
            rot_i(i) = sin1(i);
            rot_q(i) = cos1(i);
        elseif (quadrant == 1)
            rot_i(i) = cos1(i);
            rot_q(i) = sin1(i);
        elseif (quadrant == 2)
            rot_i(i) = -cos1(i);
            rot_q(i) = sin1(i);
        elseif (quadrant == 3)
            rot_i(i) = -sin1(i);
            rot_q(i) = cos1(i);
        elseif (quadrant == 4)
            rot_i(i) = sin1(i);
            rot_q(i) = -cos1(i);
        elseif (quadrant == 5)
            rot_i(i) = cos1(i);
            rot_q(i) = -sin1(i);
        elseif (quadrant == 6)
            rot_i(i) = -cos1(i);
            rot_q(i) = -sin1(i);
        elseif (quadrant == 7)
            rot_i(i) = -sin1(i);
            rot_q(i) = -cos1(i);
    end
end

dds_cos = rot_i.';
dds_sin = rot_q.';

%     spectrumScope = spectrumAnalyzer(SampleRate=20000000, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([dds_cos]);
%%
sintablen = 3200;
SINTAB = (sin(2*pi*(0:sintablen-1)./sintablen));
COSTAB = (cos(2*pi*(0:sintablen-1)./sintablen));

fs1 = 20000000;
F_required = 100000;

index = 1; 
step1 = (F_required/fs1)*sintablen;
step = round(abs(257)/16);

for i = 1:n
    dds_cos1(i) = COSTAB(round(index));
    dds_sin1(i) = SINTAB(round(index));
    index = index+step;
    if index>sintablen
        index = index-sintablen;
    end
    indexiii(i) = index;
end
figure(12)
nn = 1:n;
nn = nn.';
plot(nn, dds_cos, nn, round(dds_cos1*2^11), nn, indexiii); 

dds_cos = round(dds_cos1*2^11).';
dds_sin = round(dds_sin1*2^11).';
% 
%     spectrumScope = spectrumAnalyzer(SampleRate=fs1, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([dds_cos.']);

%%

% error_i = a(2:end)- rot_i(1:255).';
% figure(3);
% nn = (1:255);
% plot(nn, a(2:end), nn, rot_i(1:255));
% plot(error_i)

% figure(2);
% plot(rot_i)
% title('Получившийся синус из таблицы')
% xlabel('Номер отсчета') 
% ylabel('Амплитуда (int)') 

% spectrumScope = spectrumAnalyzer(SampleRate=fs, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
% spectrumScope([rot_i.']);
%%

% sintablen = 2048;
% SINTAB = sin(2*pi*(0:sintablen-1)./sintablen);
% fs1 = 2048;
% % the sintable consists of one cycle of sine wave with 2048 samples
% % if you access 2048 samples/sec (fs) from the above sintable,
% % it will generate one complete cycle of sine wave in one sec.
% % so, effectively the frequency is 1Hz.
% % step = 1; Feff(Effective frequency) = fs/sintablen;
% 
% F_required = 100;
% index = 1; step = (F_required/fs1)*sintablen;
% for i = 1:2048
%     sin1Hz(i) = SINTAB(round(index));
%     index = index+step;
%     if index>sintablen
%         index = index-sintablen;
%     end
% end
% plot(sin1Hz);
%  spectrumScope([dds_cos.']);
%%
% signal = 2*exp(j*2*pi*500000*time_base);
% % spectrumScope([signal.']);
% 
% 
% offset_sig = signal.*phase_rotation;
% % spectrumScope([offset_sig.']);
% 
% % figure(3)
% % subplot(3,1,1)
% % plot(real(phase_rotation))
% % subplot(3,1,2)
% % plot(imag(phase_rotation))
% ww = (-ang_rad/(2*16*pi/20000000)*2*pi)/20000000;
% 
% corrsignal_cordic = cordicrotate(-ang_rad,offset_sig,52);
% 
% niters_rotate = 8;
% gain = prod(sqrt(1+2.^(-2*(0:(niters_rotate-1)))));
% 
% [p_cos, p_sin] = cordic_rotate(ang_rad, gain, niters_rotate);
% 
% a = cos(ang_rad/16*time_base_1)+1i*(sin(ang_rad/16*time_base_1));
% spectrumScope([a.']);


% figure(3)
% plot(cos(p_cos*time_base_1))
% spectrumScope([corrsignal_cordic.']);
% figure(3)
% plot(imag(corrsignal_cordic))
% corrsignal = exp(-j*(ww)*time_base_1);


% out_signal1 = offset_sig.*corrsignal; %exp(-j*(ang_rad));
% x = (1:2048).';
% 
% figure(3)
% plot(x, imag(corrsignal_cordic), x, imag(out_signal1))
% 
% spectrumScope([out_signal1.']);

% nn = 14;
% x = (1:length(corrsignal_cordic)-(nn-1)).';
% z = imag(corrsignal_cordic(nn:end));
% v = imag(out_signal1(1:end-(nn-1)));
% figure(3)
% plot(x, z, x, v)
% 
% out_error1 = out_signal1(1:end-(nn-1)) - corrsignal_cordic(nn:end); 
% figure(4);
% plot(real(out_error1))








