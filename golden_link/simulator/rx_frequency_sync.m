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

    [dds_cos1, dds_sin1] = dds_int(radians_, -radians_shift, n1); 

    [dds_cos, dds_sin] = dds_int(-radians_per_sample*16, radians_per_sample, n1); 
  
    %% mult rx_signal on dds
%     out_signal_int_i = (int32(dds_cos) .* int32(real(rxsignal_int)) - int32(dds_sin) .* int32(imag(rxsignal_int)));
%     out_signal_int_q = (dds_cos .* imag(rxsignal_int) + dds_sin .* real(rxsignal_int));

%% cordic rotate
%     i_int = real(rxsignal_int);
%     q_int = imag(rxsignal_int);
%     niters_rotate = 14;
%     gain = floor(256 * prod(sqrt(1+2.^(-2*(0:(niters_rotate-1))))));
%     [p_cos, p_sin] = cordic_rotate_int(-phase_deg, gain, niters_rotate);
% 
%     i_correct = (i_int .* p_cos - q_int .* p_sin);
%     q_correct = (i_int .* p_sin + q_int .* p_cos);
% 
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
% N = length(rxsignal)-1;
% fd = 100000;
% fs = 20000000;
% bp = 12;
% bd = 16;
% ps = 90;
% % calculate frequency tuning word (FTW) from desired frequency
% ftw=round(fd/fs*2^bp);
% 
% % calculate actual synthesized frequency (Hz)
% f0=ftw*fs/2^bp;
% 
% % calculate phase tuning word (PTW) for start phase
% ptw=round(mod(ps/360*2^bp,2^bp));
% 
% % calculate actual synthesized phase start (scaled rad)
% p0=ptw/2^bp;
% 
% % create time vector for sinusoid phase (s)
% t=1/fs*[0:N-1];
% 
% % phase of sinusoid (rad)
% phase=round(mod(f0*t+p0,1)*2^bp)/2^bp*2*pi;
% 
% % quantized sinusoid value
% s=round(sin(phase)*(2^bd-1));
% 
% % normailzed sinusoid
% s=s/(2^bd-1);
%%
% n1 = 1100;
nn = (1:length(rxsignal)-1).';
% nn = (1:1100-1).';
figure(2)
a1 = dds_cos1(1:end-1);
a1 = a1;
% a11 = time_base(1:end-1).';
% a111 = a1.*a11;
a2 = real(correction_signal(2:end)*2^11);
a2 = a2.';

rx_w = importdata('test_signals\dds_out.txt');
rx_w = rx_w*2^-15;

plot(nn, a1, nn, a2);
title('Сравнение косинусов из DDS и из примера Матлаб')
xlabel('Номер отсчета') 
ylabel('Амплитуда (int)') 

error_correction = a2 - a1;

    spectrumScope = spectrumAnalyzer(SampleRate=20000000, ...            
            AveragingMethod='exponential',ForgettingFactor=0.99, ...
            YLimits=[-30 10],ShowLegend=true);

    spectrumScope([a1, a2]);
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
