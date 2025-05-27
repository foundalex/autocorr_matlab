
% Frequency error estimation and correction
function [out_signal, freq_est] = rx_frequency_sync(rxsignal, sim_options)

global sim_consts;

[n_tx_antennas, n_rx_antennas] = get_n_antennas(sim_options);

%% write to file signal
generation = 0;
if generation == 1
    ii = round(real(rxsignal)*2^11).';
    qq = round(imag(rxsignal)*2^11).';
    your_variable = [ii qq];
    dlmwrite('test_signals/rx_signal_w.txt', your_variable);
end
%%

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

    %% write signal for verilog model
    % rx_signal = rx_frequency_sync(rx_signal, sim_options);
    % i_rx = real(rx_signal);
    % q_rx = imag(rx_signal);
    % rx_signal_w = [floor(i_rx*2^11), floor(q_rx*2^11)];
    % rx_signal_w = [floor(i_rx*2^11).', floor(q_rx*2^11).'];
    % dlmwrite('test_signals/rx_signal_w.txt',rx_signal_w);
   
   %% user defined interval (last 5 symbols)
   phase1 = rxsignal(:,65:145).*conj(rxsignal(:,81:161));
%    phase1 = rxsignal(:,30:142).*conj(rxsignal(:,46:158));
   phase1 = sum(phase1, 2);
   freq_est_double = -angle(phase1) / (2*D*pi/sim_consts.SampFreq);
   angle_m = angle(phase)*180/pi;

   angle_m1 = -angle(phase1);

   %% translate to int last 5 symbols of preamb
%     rxsignal_preamb = rxsignal(2:161); 
%     rxsignal_int = int16(round(rxsignal_preamb .* 2^11).');
%     rxsignal_int_conj_q = -imag(rxsignal_int(80:160));
%     rxsignal_int_conj_i = real(rxsignal_int(80:160));
%     rxsignal_int_rt = rxsignal_int(64:144);

    width = 11;
    rxsignal_int = int16(round(rxsignal .* 2^width).');

    rxsignal_int_conj_q = -imag(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));
    rxsignal_int_conj_i = real(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));

    rxsignal_int_rt = rxsignal_int(pkt_det_offset:pkt_det_offset+rlen-D);

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
    niters_angle = 25;

    ang_cordic = cordicangle(phase,niters_angle);            % matlab function
%     ang_cordic1 = cordicangle(sum_complex,niters_angle);   % matlab function
%     ang_cordic2 = cordic_angle(phase,niters_angle);
    phase_rad = cordic_angle_int(sum_i, sum_q, niters_angle, width*2);
    er_double_phase = double(phase_rad) * 2^(-width*2); 
%     error_cordic_angle = abs(ang_cordic) - abs(er_double_phase);
%     error_cordic_angle1 = ang_cordic - ang_cordic/16;
%     error_cordic_angle2 = ang_cordic1 - ang_cordic;
    %% DDS
    n1 = length(rxsignal_int);
    phase_rad_div = bitshift(phase_rad,-17,'int32'); % divide 16
    bit_rad = dec2bin(phase_rad_div,16);
    if (bit_rad(16))
        phase_rad_div = (phase_rad_div) + 1;
    else
        phase_rad_div = (phase_rad_div);
    end
%     phase_div = round(er_double_phase*512/16);

    [dds_cos, dds_sin] = dds_int(phase_rad_div, n1);

    i_correct = int32(real(rxsignal_int)) .* int32(dds_cos) - int32(imag(rxsignal_int)) .* int32(dds_sin);
    q_correct = int32(real(rxsignal_int)) .* int32(dds_sin) + int32(imag(rxsignal_int)) .* int32(dds_cos);

    dds_cos1 = (double(dds_cos)*2^-11).';
    dds_sin1 = (double(dds_sin)*2^-11).';

    rot_mult_i1 = importdata("test_signals\rotate_mult_i.txt");
    rot_mult_i = rot_mult_i1*2^-11;
%%
figure(3)
time_base=0:n1-1;
correction_signal=exp(-j*(radians_per_sample)*time_base);
nn = (1:n1).';
a1 = rot_mult_i(1:n1);
a2 = real(correction_signal);
a3 = dds_cos1;
plot(nn, a2, nn, a3, nn, a1);

%     spectrumScope = spectrumAnalyzer(SampleRate=20000000, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([a3.']);

% subplot(2,1,1)
% plot(a1)
% title('Косинус из Verilog. Частота 100кГц')
% subplot(2,1,2)
% plot(a2)
% title('Косинуса из Матлаб. Частота 100кГц')
% xlabel('Номер отсчета') 
% ylabel('Амплитуда') 
% legend({'verilog','matlab'},'Location','southwest')

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

 out_signal = rxsignal.*correction_signal;
 i_correct_double = double(i_correct) * 2^(-width*2);
 q_correct_double = double(q_correct) * 2^(-width*2);

 out_signal1 = complex(i_correct_double, q_correct_double).';

 error_out_signal_real = abs(real(out_signal1)) - abs(real(out_signal));
 error_out_signal_imag = abs(imag(out_signal1)) - abs(imag(out_signal));
 figure(4)
 subplot(2,1,1)
 plot(error_out_signal_real)
 subplot(2,1,2)
 plot(error_out_signal_imag)
 title('Ошибка между double и int')

