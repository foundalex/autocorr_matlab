    clear all;
%%
    pkt_det_offset = 30;
    rlen = 128;
    D = 16;
    fs = 20000000;

%% Basic WLAN Link Modeling
% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 4096;          % APEP length in bytes
cfgVHT.MCS = 1;                    % Single spatial stream, 16-QAM
cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
Rs = wlanSampleRate(cfgVHT);       % Sampling rate

lstf = wlanLSTF(cfgVHT);  
lltf = wlanLLTF(cfgVHT);  
lsig = wlanLSIG(cfgVHT);
nonHTfield = [lstf;lltf;lsig]; % Combine the non-HT preamble fields
vhtsiga = wlanVHTSIGA(cfgVHT);
vhtstf = wlanVHTSTF(cfgVHT);
vhtltf = wlanVHTLTF(cfgVHT);
vhtsigb = wlanVHTSIGB(cfgVHT);
preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];
rng('default') % Initialize the random number generator
txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
data = wlanVHTData(txPSDU,cfgVHT);
txWaveform = [preamble;data]; % Transmit VHT PPDU
offset = 22.5e3;
y = frequencyOffset(txWaveform,20000000,offset);
y1 = round(y*2^11);

%% импорт данных
    rx_w = importdata('test_signals\samples.csv');
    rx_i(:,1) = rx_w(3000000:2:4000000,1);
    rx_q(:,1) = rx_w(3000000:2:4000000,2);
    rxsignal_int1 = complex(rx_i, rx_q);
%%
    rxsignal_int1 = y1;
%% ищем начало пакета
    short_preamble_detect = correlator_int(rxsignal_int1);
    index_short_preamble = find(short_preamble_detect);
    idx = index_short_preamble(1) - 99;
    rxsignal_int = rxsignal_int1(idx:end);
%% переводим в double
    rxsignal = rxsignal_int*2^-11;
%   rxsignal = rx_find_packet_edge(rxsignal2.');
%% matlab функция для поиска пакета
%     [coarsePktOffset, M] = wlanPacketDetect(rxsignal,"CBW20",0,1);
    A = wlanPacketDetect(rxsignal,"CBW20");
%% example from matlab book
    phase = rxsignal(pkt_det_offset:pkt_det_offset+rlen-D).* ...
      conj(rxsignal(pkt_det_offset+D:pkt_det_offset+rlen));
    % add all estimates 
    phase = sum(phase, 2);
    % with rx diversity combine antennas
    phase = sum(phase, 1);
    freq_est = -angle(phase) / (2*D*pi/fs);
    radians_per_sample = 2*pi*freq_est/fs;
%% example first 9 symbols
    phase1 = rxsignal(1:128).* ...
      conj(rxsignal(17:144));
    phase21 = phase1.'; 
    % add all estimates 
    phase21 = sum(phase21, 2);
    % with rx diversity combine antennas
    phase21 = sum(phase21, 1);
    freq_est21 = -angle(phase21) / (2*D*pi/fs);
%% example last 3 symbols
    phase2 = rxsignal(112:128).* ...
      conj(rxsignal(128:144));
    phase31 = phase2.'; 
    % add all estimates 
    phase31 = sum(phase31, 2);
    % with rx diversity combine antennas
    phase31 = sum(phase31, 1);
    freq_est31 = -angle(phase31) / (2*D*pi/fs);
%% matlab функция для оценки CFO
    r_ind = rxsignal(1:159);
    freqOffsetEst = wlanCoarseCFOEstimate(r_ind,'CBW20');
%% 
    % 1. Комплексно-сопряженный сигнал, задержанный на D=16 семплов
    rxsignal_int_conj_q = -imag(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));
    rxsignal_int_conj_i = real(rxsignal_int(pkt_det_offset+D:pkt_det_offset+rlen));

    rxsignal_int_rt = rxsignal_int(pkt_det_offset:pkt_det_offset+rlen-D);
%%
    % 2. Умножение задержанного комплексно-сопряженного сигнала на сигнал
    % (фазовый детектор)
    phase_detect_i = ((real(rxsignal_int_rt)) .* (rxsignal_int_conj_i) - (imag(rxsignal_int_rt)) .* (rxsignal_int_conj_q));
    phase_detect_q = ((rxsignal_int_conj_i) .* (imag(rxsignal_int_rt)) + (real(rxsignal_int_rt)) .* (rxsignal_int_conj_q));

%% 
    % sum phase
    % 3. Суммирование всех получившихся отсчетов (аккумулятор)
    sum_i = 0;
    sum_q = 0;
    for i = 1:length(rxsignal_int_rt)
        sum_i = ((sum_i) + (phase_detect_i(i)));
        sum_q = ((sum_q) + (phase_detect_q(i)));
    end
%% 
    % cordic angle
    % 4. Поиск угла получившегося значения
    niters_angle = 14;
    phase_rad = cordic_angle_int(sum_i, sum_q, niters_angle);

%%
    phase_rad_div = bitshift(phase_rad,-4,'int16'); % деление на 16
    ff = -double(phase_rad_div)*16/512/(2*pi*16*1/20000000);

    % 5. Подача на вход ДДС средне-арифметического значения (приращение фазы(шаг счетчика аккумулятора фазы))
    [dds_cos, dds_sin] = dds_int(-phase_rad_div, length(rxsignal_int));

    aa = complex(dds_cos,-dds_sin);
    aa = aa * 2^-11;

    siglen=length(rxsignal_int);
    time_base=0:siglen-1;
    correction_signal=exp(-j*(radians_per_sample)*time_base);

    spectrumScope = spectrumAnalyzer(SampleRate=fs, ...            
            AveragingMethod='exponential',ForgettingFactor=0.99, ...
            YLimits=[-30 10],ShowLegend=true);

    spectrumScope([aa.', correction_signal.']);