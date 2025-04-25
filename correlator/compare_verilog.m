
% Сравнение дефолтной модели верилога и нашей модели верилога

% Входные данные
data_i_default = importdata('test_signals\verilog\dout_i_default.txt');
data_q_default = importdata('test_signals\verilog\dout_q_default.txt');

data_i = importdata('test_signals\verilog\dout_i.txt');
data_q = importdata('test_signals\verilog\dout_q.txt');

er_i = data_i_default - data_i(1:length(data_i_default));
er_q = data_q_default - data_q(1:length(data_q_default));

%%
% Сверяем с verilog fifo
fifo_verilog = importdata('test_signals\verilog\fifo_dout.txt');
size_fifo_dat = size(fifo_verilog);
error_ver_fifo = real(delay_sample_fifo_out(18:size_fifo_dat(1)+17)) - fifo_verilog(1:end);
plot(error_ver_fifo)
%%
% Сверяем с verilog комплексный умножитель числителя
complex_mult_numerator = importdata('test_signals\verilog\complex_mult_numerator.txt');
size_ref = length(complex_mult_numerator);
ref = real(delay_prod_inst(18:size_ref+17));
error_ver_complex_numerator = ref - complex_mult_numerator(1:end);
subplot(6,1,1);
plot(error_ver_complex_numerator)
title('Ошибка после комплексного умножителя в числителе')

%% Проверка 
% сверяем с verilog окно числителя
avg_dout_dual = importdata('test_signals\verilog\avg_dout_dual.txt');
size_ref = size(avg_dout_dual);
ref = real(delay_prod_inst(18:size_ref(1)));
error_ver_avg_dual = ref - avg_dout_dual(1:size_ref(1)-17);
subplot(6,1,2);
plot(error_ver_avg_dual)
title('Ошибка после окна в числителе')
%% Проверка
% сверяем с verilog амплитуду числителя
prod_mag = importdata('test_signals\verilog\prod_mag.txt');

size_ref = size(prod_mag);
ref = mag(18:size_ref(1));

error_ver_prod_mag = ref - prod_mag(1:size_ref(1)-17);
subplot(6,1,3);
plot(error_ver_prod_mag)
title('Ошибка после вычисления амплитуды в числителе')
%% Проверка
% сверяем с verilog окно знаменателя
avg_dout = importdata('test_signals\verilog\avg_dout.txt');
size_ref = size(avg_dout);
ref = complex_mult(1:size_ref(1));
error_ver = ref - avg_dout;
subplot(6,1,5);
plot(error_ver)
title('Ошибка после окна в знаменателе')
%%  Проверка
% сверяем с verilog получившееся значение знаменателя
prod_thres  = importdata('test_signals\verilog\prod_thres.txt');

size_ref = size(prod_thres);
ref = complex_mult(18:size_ref(1));

error_ver = ref - prod_thres(1:size_ref(1)-17);
subplot(6,1,6);
plot(error_ver)
title('Ошибка в знаменателе')
%% Данные после умножителя в знаменателе 
complex_mult_mag_default = importdata('test_signals\verilog\complex_mult_mag_dout_default.txt');
complex_mult_mag = importdata('test_signals\verilog\complex_mult_mag_dout.txt');
er = complex_mult_mag_default - complex_mult_mag(1:length(complex_mult_mag_default));

% Конечные данные в знаменателе
prod_mag_default = importdata('test_signals\verilog\prod_mag_default.txt');
prod_mag = importdata('test_signals\verilog\prod_mag.txt');
er_mag = prod_mag_default - prod_mag(1:length(prod_mag_default));

% Конечные данные в числителе
avg_default = importdata('test_signals\verilog\avg_dout_default.txt');
avg = importdata('test_signals\verilog\avg_dout.txt');
er_avg = avg_default - avg(1:length(avg_default));