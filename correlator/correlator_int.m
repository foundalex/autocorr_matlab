
clear; clc

%% LOAD DATA

% wlan waveform generator
load('test_signals\var.mat');
rxWaveform = [var.waveform];

rx_w = importdata('test_signals\samples.csv');
% rx_w_24_mbps = importdata('test_signals\samples_24_mbps_A.csv');

% rx_w1(:,1) = rx_w_24_mbps(1:2:end,1);
% rx_w1(:,2) = rx_w_24_mbps(1:2:end,2);

% writematrix(rx_w1, 'test_signals\verilog\samples_blade_rf.txt');

% середина пакета
rx_blade = complex(rx_w(1:2:20000,1), rx_w(1:2:20000,2));

%% wlanPacketDetect
offset = 0;
threshold = 1;

[startOffset,M] = wlanPacketDetect(rx_blade, var.config.waveform.ChannelBandwidth, offset, threshold);
totalOffset = offset + startOffset;
plot(M)
% title(M,'')
xlabel('Отсчеты');
ylabel('Значение корреляции');
%% Генерация тестовых входных файлов для моделирования verilog
generation = 0;

switch generation
    case 1

    % нормируем сигнал
    rxWaveform_norm = rxWaveform ./ max(rxWaveform);
    
    waveform_i = floor((real(rxWaveform_norm))*2^11);
    waveform_q = floor((imag(rxWaveform_norm))*2^11);

    writematrix(waveform_i, 'test_signals\verilog\waveform_int_i.txt');
    writematrix(waveform_q, 'test_signals\verilog\waveform_int_q.txt');

    otherwise
end
%% Числитель
% Fifo
rxWaveform_norm_complex = rx_blade;
% 
% int_i = int16(real(rx_blade));
% int_q = int16(imag(rx_blade));
% 
% delay_sample_fifo_out_i = [zeros(17,1); int_i(2:end-16)];
% delay_sample_fifo_out_q = [zeros(17,1); int_q(2:end-16)];

delay_sample_fifo_out = [zeros(17,1); rxWaveform_norm_complex(2:end-16)];
% Комплексный умножитель
delay_prod_inst = rxWaveform_norm_complex .* conj(delay_sample_fifo_out);

% % conj
% delay_prod_inst_int_real = int32(int_i).* int32(delay_sample_fifo_out_i) - int32(int_q).*int32(-delay_sample_fifo_out_q);
% delay_prod_inst_int_imag = int32(int_i).* int32(delay_sample_fifo_out_q) + int32(int_i).*int32(-delay_sample_fifo_out_q);


% Отбрасываем один бит с округлением
delay_prod_inst = floor(delay_prod_inst ./ 2);
%% Проверка
% Сверяем с verilog fifo
fifo_verilog = importdata('test_signals\verilog\fifo_dout.txt');
size_fifo_dat = size(fifo_verilog);
error_ver_fifo = real(delay_sample_fifo_out(18:size_fifo_dat(1)+17)) - fifo_verilog(1:end);
plot(error_ver_fifo)

% Сверяем с verilog комплексный умножитель числителя
complex_mult_numerator = importdata('test_signals\verilog\complex_mult_numerator.txt');
size_ref = length(complex_mult_numerator);
ref = real(delay_prod_inst(18:size_ref+17));
error_ver_complex_numerator = ref - complex_mult_numerator(1:end);
subplot(6,1,1);
plot(error_ver_complex_numerator)
title('Ошибка после комплексного умножителя в числителе')
%% Окно с суммой I и Q
delay_prod = delay_prod_inst;

for i = 2 : size(delay_prod_inst)
    if (i < 34)
        delay_prod_inst(i) = delay_prod_inst(i) + delay_prod_inst(i-1);
    else
        delay_prod_inst(i) = delay_prod_inst(i-1) + delay_prod(i) - delay_prod(i-16);
    end
end

delay_prod_inst = floor(delay_prod_inst./16);
%% Проверка 
% сверяем с verilog окно числителя
avg_dout_dual = importdata('test_signals\verilog\avg_dout_dual.txt');
size_ref = size(avg_dout_dual);
ref = real(delay_prod_inst(18:size_ref(1)));
error_ver_avg_dual = ref - avg_dout_dual(1:size_ref(1)-17);
subplot(6,1,2);
plot(error_ver_avg_dual)
title('Ошибка после окна в числителе')

%% 
% abs
i_delay = abs(real(delay_prod_inst));
q_delay = abs(imag(delay_prod_inst));

% find max
for i = 1:size(delay_prod_inst)
    if (i_delay(i) > q_delay(i))
        max_d(i) = i_delay(i);
    else
        max_d(i) = q_delay(i);
    end
end

% find min
for i = 1:size(delay_prod_inst)
    if (i_delay(i) > q_delay(i))
        min_d(i) = q_delay(i);
    else
        min_d(i) = i_delay(i);
    end
end

% mag
alpha = 1;
beta = 1/4;

mag = alpha*max_d + floor(beta*min_d);
mag = mag.';

%% Проверка
% сверяем с verilog амплитуду числителя
prod_mag = importdata('test_signals\verilog\prod_mag.txt');

size_ref = size(prod_mag);
ref = mag(18:size_ref(1));

error_ver_prod_mag = ref - prod_mag(1:size_ref(1)-17);
subplot(6,1,3);
plot(error_ver_prod_mag)
title('Ошибка после вычисления амплитуды в числителе')
%% Знаменатель
% Комплексный умножитель
complex_mult = rxWaveform_norm_complex .* conj(rxWaveform_norm_complex);
% Отбрасываем один бит с округлением
complex_mult = floor(complex_mult ./ 2);

%% Проверка
% Сверяем с verilog комплексный умножитель знаменателя
complex_mult_denominator = importdata('test_signals\verilog\complex_mult_mag_dout.txt');

size_ref = size(complex_mult_denominator);
ref = complex_mult(1:size_ref(1));

error_ver_mult = ref - complex_mult_denominator;
subplot(6,1,4);
plot(error_ver_mult)
title('Ошибка после комплексного умножителя в знаменателе')
%% Окно с суммой I
% avg_channel
complex_mult_avg = complex_mult;

for i = 2 : size(complex_mult)
    if (i < 17)
        complex_mult(i) = complex_mult(i) + complex_mult(i-1);
    else
        complex_mult(i) = complex_mult(i-1) + complex_mult_avg(i) - complex_mult_avg(i-16);
    end
end

complex_mult = floor(complex_mult./16);

%% Проверка
% сверяем с verilog окно знаменателя
avg_dout = importdata('test_signals\verilog\avg_dout.txt');
size_ref = size(avg_dout);
ref = complex_mult(1:size_ref(1));
error_ver = ref - avg_dout;
subplot(6,1,5);
plot(error_ver)
title('Ошибка после окна в знаменателе')
%% 
threshold_scale = 1;
if (threshold_scale == 1)
    complex_mult = floor(complex_mult ./4) + floor(complex_mult ./8);
else
    complex_mult = floor(complex_mult ./2) + floor(complex_mult ./4);
end
%%  Проверка
% сверяем с verilog получившееся значение знаменателя
prod_thres  = importdata('test_signals\verilog\prod_thres.txt');

size_ref = size(prod_thres);
ref = complex_mult(18:size_ref(1));

error_ver = ref - prod_thres(1:size_ref(1)-17);
subplot(6,1,6);
plot(error_ver)
title('Ошибка в знаменателе')
%%
pos_count = 0;
neg_count = 0;
count = 0;

complex_mult = complex_mult(18:end);
mag = mag(18:end);

for i = 1:size(complex_mult)-(19+18)
    if (mag(i) > complex_mult(i))
        if (real(rx_blade(i+19)) < 0)
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
preamble  = importdata('test_signals\verilog\preamble.txt');
err_pream = model_pream(1:length(preamble)).' - preamble;
subplot(4,1,1);
plot(real(rx_blade(19:length(preamble))))
title('I-составляющая сигнала')
subplot(4,1,2)
plot(preamble)
title('График корреляции verilog')
subplot(4,1,3)
plot(model_pream(1:length(preamble)))
title('График корреляции matlab')
subplot(4,1,4)
plot(err_pream)
title('Ошибка корреляции matlab и verilog')