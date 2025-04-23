
% Сравнение дефолтной модели верилога и нашей модели верилога

% Входные данные
data_i_default = importdata('test_signals\verilog\dout_i_default.txt');
data_q_default = importdata('test_signals\verilog\dout_q_default.txt');

data_i = importdata('test_signals\verilog\dout_i.txt');
data_q = importdata('test_signals\verilog\dout_q.txt');

er_i = data_i_default - data_i(1:length(data_i_default));
er_q = data_q_default - data_q(1:length(data_q_default));

% Данные после умножителя в знаменателе 
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