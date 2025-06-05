function model_pream = correlator_int(rx_signal)
%% LOAD DATA
% wlan waveform generator
% load('test_signals\var.mat');
% rxWaveform = [var.waveform];
% 
% rx_w = importdata('test_signals\samples.csv');
% rx_w_24_mbps = importdata('test_signals\samples_24_mbps_A.csv');

% rx_w1(:,1) = rx_w_24_mbps(1:2:end,1);
% rx_w1(:,2) = rx_w_24_mbps(1:2:end,2);

% writematrix(rx_w1, 'test_signals\verilog\samples_blade_rf.txt');

% середина пакета
% rx_blade = complex(rx_w(1:2:20000,1), rx_w(1:2:20000,2));
% rx_blade = rx_signal;

%% wlanPacketDetect
% offset = 0;
% threshold = 1;
% 
% [startOffset,M] = wlanPacketDetect(rx_blade, var.config.waveform.ChannelBandwidth, offset, threshold);
% totalOffset = offset + startOffset;
% plot(M)
% % title(M,'')
% xlabel('Отсчеты');
% ylabel('Значение корреляции');

%% Числитель
% Fifo
int_i = real(rx_signal);
int_q = imag(rx_signal);

delay_sample_fifo_out_i = [zeros(17,1); int_i(2:end-16)];
delay_sample_fifo_out_q = [zeros(17,1); int_q(2:end-16)];

delay_sample_fifo_out = [zeros(17,1); rx_signal(2:end-16)];

% Комплексный умножитель
rxsignal_int_conj_q = -imag(delay_sample_fifo_out);
rxsignal_int_conj_i = real(delay_sample_fifo_out);

% Комплексное умножение на комплексно-сопряженное число
delay_prod_inst_i = (int32(real(rx_signal)) .* int32(rxsignal_int_conj_i) - int32(imag(rx_signal)) .* int32(rxsignal_int_conj_q));
delay_prod_inst_q = (int32(rxsignal_int_conj_i) .* int32(imag(rx_signal)) + int32(real(rx_signal)) .* int32(rxsignal_int_conj_q));

% Отбрасываем один бит с округлением
delay_prod_inst_i = bitshift(delay_prod_inst_i,-1,'int32');
delay_prod_inst_q = bitshift(delay_prod_inst_q,-1,'int32');

% Сверяем результаты
delay_prod_inst_i_v = importdata('test_signals\sync_short\delay_prod_i.txt');
delay_prod_inst_q_v = importdata('test_signals\sync_short\delay_prod_q.txt');
error_delay_prod_inst_i = delay_prod_inst_i(18:length(delay_prod_inst_i_v)+17) - int32(delay_prod_inst_i_v);
error_delay_prod_inst_q = delay_prod_inst_q(18:length(delay_prod_inst_q_v)+17) - int32(delay_prod_inst_q_v);

%% Окно с суммой I и Q
delay_prod_i = delay_prod_inst_i;
delay_prod_q = delay_prod_inst_q;

for i = 2 : size(delay_prod_inst_i)
    if (i < 34)
        delay_prod_inst_i(i) = delay_prod_inst_i(i) + delay_prod_inst_i(i-1);
        delay_prod_inst_q(i) = delay_prod_inst_q(i) + delay_prod_inst_q(i-1);
    else
        delay_prod_inst_i(i) = delay_prod_inst_i(i-1) + delay_prod_i(i) - delay_prod_i(i-16);
        delay_prod_inst_q(i) = delay_prod_inst_q(i-1) + delay_prod_q(i) - delay_prod_q(i-16);
    end
end

% делим на 16
delay_prod_inst_i = bitshift(delay_prod_inst_i,-4,'int32');
delay_prod_inst_q = bitshift(delay_prod_inst_q,-4,'int32');

% Сверяем результаты
prod_avg_i_v = importdata('test_signals\sync_short\prod_avg_i.txt');
prod_avg_q_v = importdata('test_signals\sync_short\prod_avg_q.txt');
error_prod_avg_i = delay_prod_inst_i(18:length(prod_avg_i_v)+17) - int32(prod_avg_i_v);
error_prod_avg_q = delay_prod_inst_q(18:length(prod_avg_q_v)+17) - int32(prod_avg_q_v);

%% 
% abs
i_delay = abs(delay_prod_inst_i);
q_delay = abs(delay_prod_inst_q);

% find max
for i = 1:length(delay_prod_inst_i)
    if (i_delay(i) > q_delay(i))
        max_d(i) = i_delay(i);
    else
        max_d(i) = q_delay(i);
    end
end

% find min
for i = 1:length(delay_prod_inst_i)
    if (i_delay(i) > q_delay(i))
        min_d(i) = q_delay(i);
    else
        min_d(i) = i_delay(i);
    end
end

% mag
alpha = 1;
beta = 1/4;

mag = max_d + bitshift(min_d,-2,'int32');
mag = mag.';

% Сверяем результаты
prod_avg_mag_v = importdata('test_signals\sync_short\prod_avg_mag.txt');
error_prod_avg_mag = mag(18:length(prod_avg_mag_v)+17) - int32(prod_avg_mag_v);

%% Знаменатель
rx_signal_int_conj_i = real(rx_signal);
rx_signal_int_conj_q = -imag(rx_signal);

complex_mult_i = (int32(real(rx_signal)) .* int32(rx_signal_int_conj_i) - ...
    int32(imag(rx_signal)) .* int32(rx_signal_int_conj_q));
complex_mult_q = (int32(rx_signal_int_conj_i) .* int32(imag(rx_signal)) + ...
    int32(real(rx_signal)) .* int32(rx_signal_int_conj_q));

complex_mult_i = bitshift(complex_mult_i,-1,'int32');
complex_mult_q = bitshift(complex_mult_q,-1,'int32');

%% Окно с суммой
% avg_channel
complex_mult_avg_i = complex_mult_i;

for i = 2 : length(complex_mult_i)
    if (i < 17)
        complex_mult_i(i) = complex_mult_i(i) + complex_mult_i(i-1);
    else
        complex_mult_i(i) = complex_mult_i(i-1) + complex_mult_avg_i(i) - complex_mult_avg_i(i-16);
    end
end

complex_mult_i = bitshift(complex_mult_i,-4,'int32');
%% 
threshold_scale = 1;
if (threshold_scale == 1)
    complex_mult_i = bitshift(complex_mult_i,-2,'int32') + bitshift(complex_mult_i,-3,'int32');
else
    complex_mult_i = bitshift(complex_mult_i,-1,'int32') + bitshift(complex_mult_i,-2,'int32');
end

%%
pos_count = 0;
neg_count = 0;
count = 0;

mag = mag(18:end);

for i = 1:length(complex_mult_i)-(19+18)
    if (mag(i) > complex_mult_i(i))
        if (int_i(i+19) < 0)
            neg_count = neg_count + 1;
        else
            pos_count = pos_count + 1;
        end
    
        if (count > 100)
            if (pos_count > 25 && neg_count > 25)
                model_pream(i) = 1;
            end
            count = 0;
            pos_count = 0;
            neg_count = 0;
        else 
            count = count + 1;
        end
    else 
        model_pream(i) = 0;
        pos_count = 0;
        neg_count = 0;
        count = 0;
    end
end
%% 
% preamble  = importdata('test_signals\verilog\preamble.txt');
% err_pream = model_pream(1:length(preamble)).' - preamble;
% subplot(4,1,1);
% plot(real(rx_blade(19:length(preamble))))
% title('I-составляющая сигнала')
% subplot(4,1,2)
% plot(preamble)
% title('График корреляции verilog')
% subplot(4,1,3)
% plot(model_pream(1:length(preamble)))
% title('График корреляции matlab')
% subplot(4,1,4)
% plot(err_pream)
% title('Ошибка корреляции matlab и verilog')