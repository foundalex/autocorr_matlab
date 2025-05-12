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
   
   % add all estimates 
   phase = sum(phase, 2);
   
   % with rx diversity combine antennas
   phase = sum(phase, 1);
   
   freq_est = -angle(phase) / (2*D*pi/sim_consts.SampFreq);
   
   radians_per_sample = 2*pi*freq_est/sim_consts.SampFreq;

   %% matlab function 
   freq_est_matlab = wlanCoarseCFOEstimate(rxsignal(2:161).','CBW20');  % if no packet offset

   %% user defined interval (last 5 symbols)
   phase1 = rxsignal(:,65:145).*conj(rxsignal(:,81:161));
%    phase1 = rxsignal(:,30:142).*conj(rxsignal(:,46:158));
   phase1 = sum(phase1, 2);
   freq_est_double = -angle(phase1) / (2*D*pi/sim_consts.SampFreq);
   ang_cordic = cordicangle(phase1,8);
   ang_cordic_user = cordic_angle(phase1,8);
   
   %% cordic angle

   sin_table = [0.707; 0.4472; 0.2425; 0.124; 0.0624; 0.0312; 0.0156];
   cos_table = [0.707; 0.8944; 0.9701; 0.9923; 0.9981; 0.9995; 0.9999];

   tan_v = imag(phase)/real(phase);
   xn = 1;
   yn = 0;

   for i = 1:6
       if (yn < 0)
            xn = xn*cos_table(i) - yn*sin_table(i);
            yn = yn*cos_table(i) + xn*sin_table(i);
       else
            xn = xn*cos_table(i) + yn*sin_table(i);
            yn = yn*cos_table(i) - xn*sin_table(i);
       end
   end

   %% int 
    
    %debug
%     phase_d = rxsignal(:,65:145).*conj(rxsignal(:,81:161));
%     phase_d = phase_d.' .* 2^22;

    %translate to int
    rxsignal_preamb = rxsignal(2:161); 
    rxsignal_int = int16(round(rxsignal_preamb .* 2^11).');
    rxsignal_int_conj_q = -imag(rxsignal_int(80:160));
    rxsignal_int_conj_i = real(rxsignal_int(80:160));

    rxsignal_int_rt = rxsignal_int(64:144);

    phase_detect_i = int32(real(rxsignal_int_rt)) .* int32(rxsignal_int_conj_i) - int32(imag(rxsignal_int_rt)) .* int32(rxsignal_int_conj_q);
    phase_detect_q = int32(rxsignal_int_conj_i) .* int32(imag(rxsignal_int_rt)) + int32(real(rxsignal_int_rt)) .* int32(rxsignal_int_conj_q);

%     phase_detect_complex = complex(phase_detect_i,phase_detect_q);
%     phase_detect_sum = sum(phase_detect_complex,1);
    
    sum_i = 0;
    sum_q = 0;
    for i = 1:length(rxsignal_int_rt)
        sum_i = sum_i + phase_detect_i(i);
        sum_q = sum_q + phase_detect_q(i);
    end

    samples = (0:length(rxsignal)-1);

    ang = cordicatan2(sum_q,sum_i);
    denominator = (2*D*pi/sim_consts.SampFreq);
    %%
%     freq_est_int = -double(ang) / denominator;
%     correction_signal_i = cos(2*pi*freq_est_int*samples/sim_consts.SampFreq);
%     correction_signal_q = sin(2*pi*freq_est_int*samples/sim_consts.SampFreq);
%     correction_signal_complex = complex(correction_signal_i,correction_signal_q);
    %%
    % from openofdm 
    freq_est_ofdm = double(ang) / 16;
    freq_ofdm = -freq_est_ofdm / (2 * pi/20000000);
%     correction_signal_i_ofdm = cos(freq_est_ofdm*samples);
%     correction_signal_q_ofdm = sin(freq_est_ofdm*samples);
%     correction_signal_complex_ofdm = complex(correction_signal_i_ofdm,correction_signal_q_ofdm);

%     freq_est_ofdm_int = ang / 16;
%     correction_signal_i_ofdm_int = cos(freq_est_ofdm*samples);
%     correction_signal_q_ofdm_int = sin(freq_est_ofdm*samples);
%     correction_signal_complex_ofdm_int = complex(correction_signal_i_ofdm_int,correction_signal_q_ofdm_int);

    freq_err = freq_est - freq_ofdm;

    out_signal_int = cordicrotate(freq_est_ofdm,rxsignal);
    %%
    % error openofdm and example 
%     err = correction_signal_complex - complex(correction_signal_i_ofdm,-correction_signal_q_ofdm);
%     figure(2);
%     plot(imag(err))

else
   % Magic number
   freq_est = sim_options.FreqError;
   radians_per_sample = 2*pi*freq_est/sim_consts.SampFreq;
end

% Now create a signal that has the frequency offset in the other direction
siglen=length(rxsignal(1,:));
time_base=0:siglen-1;
correction_signal=repmat(exp(-j*(radians_per_sample)*time_base),n_rx_antennas,1);

% % compare two methods
% error = correction_signal - correction_signal_complex; 
% figure(3);
% plot(real(error))

% And finally apply correction on the signal

out_signal = rxsignal.*correction_signal;

out_error = out_signal - out_signal_int; 
figure(2);
subplot(2,1,1)
plot(real(out_error))
subplot(2,1,2)
plot(imag(out_error))

