% Frequency error estimation and correction
function [out_signal, freq_est] = rx_frequency_sync(rxsignal, sim_options)

global sim_consts;

[n_tx_antennas, n_rx_antennas] = get_n_antennas(sim_options);

%% write to file signal
ii = round(real(rxsignal)*2^11).';
qq = round(imag(rxsignal)*2^11).';
your_variable = [ii qq];
dlmwrite('test_signals/rx_signal_w.txt', your_variable);
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

   %% matlab function 
%    freq_est_matlab = wlanCoarseCFOEstimate(rxsignal(2:161).','CBW20');  % if no packet offset

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
%    ang_cordic_user = cordic_angle(phase1,8);

   %% translate to int last 5 symbols of preamb
%     rxsignal_preamb = rxsignal(2:161); 
%     rxsignal_int = int16(round(rxsignal_preamb .* 2^11).');
%     rxsignal_int_conj_q = -imag(rxsignal_int(80:160));
%     rxsignal_int_conj_i = real(rxsignal_int(80:160));
%     rxsignal_int_rt = rxsignal_int(64:144);

    rxsignal_int = int16(round(rxsignal .* 2^11).');

    rxsignal_int_conj_q = -imag(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));
    rxsignal_int_conj_i = real(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));

    rxsignal_int_rt = rxsignal_int(pkt_det_offset:pkt_det_offset+rlen-D);

    phase_detect_i = int32(real(rxsignal_int_rt)) .* int32(rxsignal_int_conj_i) - int32(imag(rxsignal_int_rt)) .* int32(rxsignal_int_conj_q);
    phase_detect_q = int32(rxsignal_int_conj_i) .* int32(imag(rxsignal_int_rt)) + int32(real(rxsignal_int_rt)) .* int32(rxsignal_int_conj_q);

    %% error phase
    phase_detect_i_double = double(phase_detect_i) * 2^-22;
    phase_detect_q_double = double(phase_detect_q) * 2^-22;

    error_ph_i = real(phase2) - phase_detect_i_double;
    error_ph_q = imag(phase2) - phase_detect_q_double;
    %% sum phase
    sum_i = 0;
    sum_q = 0;
    for i = 1:length(rxsignal_int_rt)
        sum_i = sum_i + phase_detect_i(i);
        sum_q = sum_q + phase_detect_q(i);
    end
    %% error sum
    sum_i_double = double(sum_i) * 2^-22;
    sum_q_double = double(sum_q) * 2^-22;
    sum_complex = complex(sum_i_double,sum_q_double);
    error_sum = phase - sum_complex;
    %% cordic angle
    niters_angle = 14;
    ang_cordic = cordicangle(phase,niters_angle); % matlab function
    ang_cordic1 = cordicangle(sum_complex,niters_angle); % matlab function
%     ang_cordic2 = cordic_angle(sum_complex,8);

    phase_rad = cordic_angle_int(sum_i, sum_q, 14);

%     ang_error = angle_m - (phase_deg/256); % compare results
%     phase_rad_shifted = (bitshift(-phase_rad,-4,'int32')); % window size 16
%     cordic_user_rad = -(double((phase_deg/256))*pi/180);
%     phase_deg_example = radians_per_sample * 180/pi;
%     a = double(phase_deg)/256*pi/180
%     error_angle = radians_per_sample - a;

    n1 = length(rxsignal_int);

    radians_ = double(phase_rad)/512;
    radians_shift = double(phase_rad)/16/512;
%     [dds_cos2, dds_sin2] = dds_int_my(phase_rad,n1);

%     [dds_cos2, dds_sin2] = dds_int(radians_, -radians_shift, n1); 
% 
%     [dds_cos, dds_sin] = dds_int(-radians_per_sample*16, radians_per_sample, n1); 
  
    %% mult rx_signal on dds
%     out_signal_int_i = (int32(dds_cos) .* int32(real(rxsignal_int)) - int32(dds_sin) .* int32(imag(rxsignal_int)));
%     out_signal_int_q = (dds_cos .* imag(rxsignal_int) + dds_sin .* real(rxsignal_int));

    %% cordic rotate
    i_int = real(rxsignal_int);
    q_int = imag(rxsignal_int);
    niters_rotate = 14;
    pi_int = round(pi*512);
    gain = floor(512 * prod(sqrt(1+2.^(-2*(0:(niters_rotate-1))))));
    sig = cordicrotate(radians_per_sample,rxsignal,14);
    phaser = 0;
    for i = 1:length(rxsignal_int)    
        [p_cos(i), p_sin(i)] = cordic_rotate_int(-phaser, gain, niters_rotate);

        phaser = phaser + round(phase_rad/16);
        if (phaser < -pi_int)
            phaser = phaser + 2*pi_int;
        elseif (phaser > pi_int)
            phaser = phaser - 2*pi_int;
        end
    end
%         i_correct = (i_int .* p_cos - q_int .* p_sin);
%         q_correct = (i_int .* p_sin + q_int .* p_cos);
%     end

%     i_correct_double = double(i_correct).*2^-11;
%     q_correct_double = double(q_correct).*2^-11;
%     corrected_signal_complex = complex(i_correct_double,q_correct_double).';

    %% OPEN OFDM
%     ang = cordicatan2(sum_q,sum_i);
%     freq_est_ofdm = double(ang) / 16;
%     freq_ofdm = -freq_est_ofdm / (2 * pi/20000000);
%     freq_err = freq_est - freq_ofdm;
% 
%     out_signal_int = cordicrotate(freq_est_ofdm,rxsignal);

else
   % Magic number
   freq_est = sim_options.FreqError;
   radians_per_sample = 2*pi*freq_est/sim_consts.SampFreq;
end

% Now create a signal that has the frequency offset in the other direction

siglen=length(rxsignal(1,:));
time_base=0:siglen-1;
correction_signal=repmat(exp(-j*(radians_per_sample)*time_base),n_rx_antennas,1);
%%
% sintablen = 2^11;
% SINTAB = (-sin(2*pi*(0:sintablen-1)./sintablen));
% COSTAB = (cos(2*pi*(0:sintablen-1)./sintablen));
% % 
% % SINTAB = round(SINTAB*sintablen);
% % COSTAB = round(COSTAB*sintablen);
% 
% fs1 = 20000000;
% F_required = 100000;
% index = 1; 
% index1 = 1;
% % step = ((F_required/fs1)*sintablen);
% step = 16;
% step1 = round((F_required/fs1)*sintablen);
% 
% for i = 1:n1
%     dds_cos1(i) = COSTAB((index));
%     dds_sin1(i) = SINTAB(round(index));
% 
%     dds_cos2(i) = COSTAB(round(index1));
%     dds_sin2(i) = SINTAB(round(index1));
% 
%     index = index+step;
%     index1 = index1+step1;
% 
%     if index>sintablen
%       index = index - sintablen;
%     end
% 
%     if index1>sintablen
%       index1 = index1 - sintablen;
%     end
%     indexiii(i) = index;
%     indexiii1(i) = index1;
% end

% figure(13)
% n2 = (1:n1).';
% plot(n2, indexiii, n2, indexiii1)

% figure(14)
% subplot(2,1,1)
% plot(COSTAB);
% subplot(2,1,2)
% plot(angle(complex(COSTAB,SINTAB)));
%%
% lines = readlines("test_signals\rot_lut.txt");
% % a = importdata("test_signals\rotate_mult_i.txt");
% 
% % convert int sine to double
% for i = 1:512
%     str = num2str(lines(i));
%     r = str - '0';
%     rr = dec2bin(r);
%     str_q = [rr(17), rr(18), rr(19), rr(20), rr(21), rr(22), rr(23), rr(24), rr(25), rr(26), rr(27), rr(28), rr(29), rr(30), rr(31), rr(32)];
%     str_i = [rr(1), rr(2), rr(3), rr(4), rr(5), rr(6), rr(7), rr(8), rr(9), rr(10), rr(11), rr(12), rr(13), rr(14), rr(15), rr(16)];
% 
%     str_i_array(i) = bin2dec(str_i);
%     str_q_array(i) = bin2dec(str_q);
% end
% qw = (1:512).';
% figure(20)
% plot(qw, str_i_array, qw, round(COSTAB(1:512)*2048))
% 
% % dds_cos1 = round(dds_cos1*2^11).';
% % dds_sin1 = round(dds_sin1*2^11).';


%%
nn = (1:length(rxsignal)).';
figure(12);
a1 = p_cos;
a1 = a1;
% a2 = real(correction_signal*2^11);
a2 = real(correction_signal);
a2 = a2.';
as = angle(correction_signal);
% ar = angle(complex(dds_cos1,dds_sin1));
% plot(nn, ar, nn, as);


rx_w = importdata('test_signals\dds_out.txt');
rx_w = rx_w*2^-15;
figure(2);
plot(nn, a1, nn, a2);
title('Сравнение косинусов из DDS и из примера Матлаб')
xlabel('Номер отсчета') 
ylabel('Амплитуда (int)') 

error_correction = a2 - a1;

%     spectrumScope = spectrumAnalyzer(SampleRate=1000000, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([a1.', a2]);
% spectrumScope([rx_w]);
%% And finally apply correction on the signal

 out_signal = rxsignal.*correction_signal;

%%
signal_int = int32(floor(real(out_signal(1:1024)*2^11))).';
out_error = signal - out_signal_int_i;



% figure(2);
% subplot(2,1,1)
% plot(real(out_error))
% title('Ошибка I-составляющей')
% subplot(2,1,2)
% plot(imag(out_error))
% title('Ошибка Q-составляющей')

%%


% figure(3)
% subplot(3,1,1)
% plot(real(phase_rotation))
% subplot(3,1,2)
% plot(imag(phase_rotation))

% corrsignal_cordic = cordicrotate(ang_rad,offset_sig,52);
% 
% % figure(3)
% % plot(imag(corrsignal_cordic))
% 
% corrsignal = exp(-j*(ang_rad)*time_base1);
% out_signal1 = offset_sig.*corrsignal;
% 
% nn = 14;
% x = (1:length(corrsignal_cordic)-(nn-1)).';
% 
% z = imag(corrsignal_cordic(nn:end));
% 
% v = imag(out_signal1(1:end-(nn-1)));
% figure(3)
% plot(x, z, x, v)
% 
% 
% out_error1 = out_signal1(1:end-(nn-1)) - corrsignal_cordic(nn:end); 
% figure(4);
% plot(real(out_error1))
% figure(4)
% plot(real(corrsignal_cordic))
% figure(5)
% plot(real(out_signal1))
