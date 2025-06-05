
% Frequency error estimation and correction
function [out_signal, freq_est] = rx_frequency_sync(rxsignal, sim_options)

global sim_consts;

[n_tx_antennas, n_rx_antennas] = get_n_antennas(sim_options);

% Estimate the frequency error
if sim_options.FreqSync
   
   % allows for error in packet detection
   pkt_det_offset = 30;
   
   % averaging length
   rlen = 128;
   
   % short training symbol periodicity
   D = 16;
   
   phase = rxsignal(:,pkt_det_offset:pkt_det_offset+rlen-D).* ...
      conj(rxsignal(:,pkt_det_offset+D:pkt_det_offset+rlen));

   phase2 = phase.'; 
   
   % add all estimates 
   phase = sum(phase, 2);
   
   % with rx diversity combine antennas
   phase = sum(phase, 1);
   
   freq_est = -angle(phase) / (2*D*pi/sim_consts.SampFreq);
   
   radians_per_sample = 2*pi*freq_est/sim_consts.SampFreq;

   %% user defined interval (last 5 symbols)
   phase1 = rxsignal(:,65:145).*conj(rxsignal(:,81:161));
%    phase1 = rxsignal(:,30:142).*conj(rxsignal(:,46:158));
   phase1 = sum(phase1, 2);
   freq_est_double = -angle(phase1) / (2*D*pi/sim_consts.SampFreq);
   angle_m = angle(phase)*180/pi;
   angle_m1 = -angle(phase1);

   width = 11;
   rxsignal_int = int16(round(rxsignal .* 2^width).');
   %% write to file signal
   generation = 0;

   if generation == 1
    ii = real(rxsignal_int);
    qq = imag(rxsignal_int);
    your_variable = [ii qq];
    writematrix(your_variable, 'test_signals/rx_signal_fo_17000.txt');

    rx_w = importdata('test_signals\samples.csv');
    rx_w1(:,1) = rx_w(1:2:length(rx_w)/10,1);
    rx_w1(:,1) = rx_w(1:2:length(rx_w)/10,2);
    writematrix(rx_w1, 'test_signals\verilog\samples_blade_rf.txt');
   end

    % Поиск короткой преамбулы
    short_preamble_detect = correlator_int(rxsignal_int);

    index_short_preamble = find(short_preamble_detect);
    rxsignal_int = rxsignal_int(index_short_preamble-101:end);

    % 1. Комплексно-сопряженный сигнал, задержанный на D=16 семплов
    rxsignal_int_conj_q = -imag(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));
    rxsignal_int_conj_i = real(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));

    rxsignal_int_rt = rxsignal_int(pkt_det_offset:pkt_det_offset+rlen-D);

    % 2. Умножение задержанного комплексно-сопряженного сигнала на сигнал
    % (фазовый детектор)
    phase_detect_i = (int32(real(rxsignal_int_rt)) .* int32(rxsignal_int_conj_i) - int32(imag(rxsignal_int_rt)) .* int32(rxsignal_int_conj_q));
    phase_detect_q = (int32(rxsignal_int_conj_i) .* int32(imag(rxsignal_int_rt)) + int32(real(rxsignal_int_rt)) .* int32(rxsignal_int_conj_q));

    phase_detect_complex_double = complex(double(rxsignal_int_conj_i),double(rxsignal_int_conj_q))*2^-width .* (double(rxsignal_int_rt))*2^-width;
    error_phase_ph = phase2 - phase_detect_complex_double;

    %% error phase
    phase_detect_i_double = double(phase_detect_i) * 2^-(width*2);
    phase_detect_q_double = double(phase_detect_q) * 2^-(width*2);

    error_ph_i = real(phase2) - phase_detect_i_double;
    error_ph_q = imag(phase2) - phase_detect_q_double;
    %% sum phase
    % 3. Суммирование всех получившихся отсчетов (аккумулятор)
    sum_i = 0;
    sum_q = 0;
    for i = 1:length(rxsignal_int_rt)
        sum_i = (int32(sum_i) + int32(phase_detect_i(i)));
        sum_q = (int32(sum_q) + int32(phase_detect_q(i)));
    end
    %% error sum
    sum_i_double = double(sum_i) * 2^(-width*2);
    sum_q_double = double(sum_q) * 2^(-width*2);
    sum_complex = complex(sum_i_double,sum_q_double);
    error_sum = phase - sum_complex;
    %% cordic angle
    % 4. Поиск угла получившегося значения
    niters_angle = 7;
    ang_cordic = cordicangle(phase,niters_angle);            % matlab function
%     ang_cordic1 = cordicangle(sum_complex,niters_angle);   % matlab function
%     ang_cordic2 = cordic_angle(phase,niters_angle);
    phase_rad = cordic_angle_int(sum_i, sum_q, niters_angle);
    er_double_phase = double(phase_rad) * 2^(-width*2); 
%     error_cordic_angle = abs(ang_cordic) - abs(er_double_phase);
%     error_cordic_angle1 = ang_cordic - ang_cordic/16;
%     error_cordic_angle2 = ang_cordic1 - ang_cordic;
    %% DDS
    % 5. Нахождение среднего арифметического (деление на 16)
    n1 = length(rxsignal_int);
    bit_rad = dec2bin(phase_rad,16);
    phase_rad_div = bitshift(phase_rad,-4,'int16'); % деление на 16
    bit_rad1 = dec2bin(phase_rad_div,16);

    % 6. Подача на вход ДДС средне-арифметического значения (приращение фазы(шаг счетчика аккумулятора фазы))
    [dds_cos, dds_sin] = dds_int(phase_rad_div, n1+2000);

    dds_cos2 = dds_cos(1:n1).';
    dds_sin2 = dds_sin(1:n1).';
    % 7. Смещение частоты входного сигнала путем умножения на выходной
    % сигнал ДДС
    i_correct = int32(real(rxsignal_int)) .* int32(dds_cos2) - int32(imag(rxsignal_int)) .* int32(dds_sin2);
    q_correct = int32(real(rxsignal_int)) .* int32(dds_sin2) + int32(imag(rxsignal_int)) .* int32(dds_cos2);

    dds_cos1 = (double(dds_cos)*2^-11).';
    dds_sin1 = (double(dds_sin)*2^-11).';

    rot_mult_cos = importdata("test_signals\rotate_mult_i.txt");
    rot_mult_i = rot_mult_cos*2^-11;
    rot_mult_sin = importdata("test_signals\rotate_mult_q.txt");
    rot_mult_q = rot_mult_sin*2^-11;

else
   % Magic number
   freq_est = sim_options.FreqError;
   radians_per_sample = 2*pi*freq_est/sim_consts.SampFreq;
end

% Now create a signal that has the frequency offset in the other direction

siglen=length(rxsignal(1,:));
time_base=0:siglen-1;
correction_signal=repmat(exp(-j*(radians_per_sample)*time_base),n_rx_antennas,1);
%% And finally apply correction on the signal

% figure(3)
% time_base=0:n1-1;
% correction_signal=exp(-j*(radians_per_sample)*time_base);
% nn = (1:n1).';
% a1 = rot_mult_i(1:n1);
% a2 = real(correction_signal);
% a3 = dds_cos1(1:n1);
% plot(nn, a2, nn, a3, nn, a1);
% 
% figure(8)
% error_correction = a2 - a3;
% plot(error_correction);
% title('Разница между генерируемыми косинусами матлаб и verilog')
% xlabel('Номер отсчета') 
% ylabel('Величина ошибки') 
% % legend({'verilog','matlab'},'Location','southwest')

%     spectrumScope = spectrumAnalyzer(SampleRate=sim_consts.SampFreq, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([dds_cos2]);

 out_signal = rxsignal.*correction_signal;

 %% Ошибка между моей моделью и матлаб
 i_correct_double = double(i_correct) * 2^((-width*2));
 q_correct_double = double(q_correct) * 2^((-width*2));

 out_signal1 = complex(i_correct_double, q_correct_double);

 error_out_signal_real = abs(real(out_signal1)) - abs(real(out_signal).');
 error_out_signal_imag = abs(imag(out_signal1)) - abs(imag(out_signal).');
 figure(4)
 subplot(2,1,1)
 plot(error_out_signal_real)
 title('Ошибка между I-составляющими сигнала после смесителя (алгоритм matlab double и мой int)')
 subplot(2,1,2)
 plot(error_out_signal_imag)
 title('Ошибка между Q-составляющими сигнала после смесителя (алгоритм matlab double и мой int)')
 %% Ошибка между верилог моделью и своей моделью
 % необходимо на этом месте остановить моделирование, чтобы сигнал не менялся из-за шума
 % *************************************
 rotate_out_i = importdata("test_signals\rotate_out_i.txt");
 rotate_out_i1 = rotate_out_i*2^-22;
 rotate_out_q = importdata("test_signals\rotate_out_q.txt");
 rotate_out_q1 = rotate_out_q*2^-22;

 count = 1;
 count1 = 1000;
 for kk = 1:length(rxsignal_int)
    aa(kk) = int32(dds_cos(kk));
    bb(kk) = int32(dds_sin(kk));
    cc(kk) = int32(real(rxsignal_int(count+189)));
    dd(kk) = int32(imag(rxsignal_int(count+189)));
    i_correct_verilog(kk) = cc(kk) .* aa(kk) - dd(kk) .* bb(kk);
    i_correct_verilog(kk) = bitshift(i_correct_verilog(kk),-1,'int32');
    q_correct_verilog(kk) = cc(kk) .* bb(kk) + dd(kk) .* aa(kk);
    q_correct_verilog(kk) = bitshift(q_correct_verilog(kk),-1,'int32');
    count = count + 1;
    if (count == 129 || count == (count1+64))
        count = count + 16;
        count1 = count;
    end

    if (count > length(rxsignal_int)-190)
        break;
    end
 end

%  i_correct_verilog_double = double(i_correct_verilog).' * 2^(-width*2);
%  q_correct_verilog_double = double(q_correct_verilog).' * 2^(-width*2);
%  error_verilog_out_i = abs(rotate_out_i1(1:2200)) - abs(i_correct_verilog_double(1:2200));
%  error_verilog_out_q = abs(rotate_out_q1(1:2200)) - abs(q_correct_verilog_double(1:2200));

 i_correct_verilog = i_correct_verilog.';
 q_correct_verilog = q_correct_verilog.'; 
 error_verilog_out_i = abs(rotate_out_i(1:2200)) ./ abs(double(i_correct_verilog(1:2200)));
 error_verilog_out_q = abs(rotate_out_q(1:2200)) ./ abs(double(q_correct_verilog(1:2200)));

 figure(5)
 subplot(2,1,1)
 plot(error_verilog_out_i)
 title('Отношение между отсчетами I-составляющей')
 xlabel('Номер отсчета') 
 ylabel('Величина ошибки') 
 subplot(2,1,2)
 plot(error_verilog_out_q)
 title('Отношение между отсчетами Q-составляющей')
 xlabel('Номер отсчета') 
 ylabel('Величина ошибки') 
%  legend({'verilog','matlab'},'Location','northeast')

    spectrumScope = spectrumAnalyzer(SampleRate=sim_consts.SampFreq, ...            
            AveragingMethod='exponential',ForgettingFactor=0.99, ...
            YLimits=[-30 10],ShowLegend=true);
    ar = complex(dds_cos2,dds_sin2)
    spectrumScope([ar]);