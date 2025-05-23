clear all;

VAL = 2*pi;


ATAN_LUT_SCALE = 512;
SCALE = 2048;

MAX1 = round(VAL*ATAN_LUT_SCALE);
SIZE = 2^ceil(log2(MAX1));
%%
lines = readlines("test_signals\rot_lut.txt");
% convert int sine to double
for i = 1:512
    str = num2str(lines(i));
    r = str - '0';
    rr = dec2bin(r);
    str_q = [rr(17), rr(18), rr(19), rr(20), rr(21), rr(22), rr(23), rr(24), rr(25), rr(26), rr(27), rr(28), rr(29), rr(30), rr(31), rr(32)];
    str_i = [rr(1), rr(2), rr(3), rr(4), rr(5), rr(6), rr(7), rr(8), rr(9), rr(10), rr(11), rr(12), rr(13), rr(14), rr(15), rr(16)];

    str_i_array(i) = bin2dec(str_i);
    str_q_array(i) = bin2dec(str_q);
end

%% 
for i = 1:SIZE
    key(i) = i/MAX1*VAL;
    ii(i) = round(cos(key(i))*SCALE);
    qq(i) = round(-sin(key(i))*SCALE);
end

i1 = int32(ii*2^16);
i2 = ii/2048;
q2 = qq/2048;

index = 1;
step = 16;
COSTAB = cos(2*pi*(0:SIZE-1)./SIZE);
SINTAB = -sin(2*pi*(0:SIZE-1)./SIZE);

figure(13)
nn = (1:SIZE).';
plot(nn, i2, nn, COSTAB);
%% verilog table
cos_table = cos(key);
sin_table = sin(key);

atan_ans1 = atan(sin_table(1)/cos_table(1));
atan_ans2 = atan(sin_table(2)/cos_table(2));

ff1 = -atan_ans1/(2*pi/200000000);
ff2 = -atan_ans2/(2*pi/200000000);

% matlab table
atan_ans3 = atan(SINTAB(2)/COSTAB(2));
atan_ans4 = atan(SINTAB(3)/COSTAB(3));

ff3 = -atan_ans3/(2*pi/200000000);
ff4 = -atan_ans4/(2*pi/200000000);

%     spectrumScope = spectrumAnalyzer(SampleRate=1000000, ...            
%             AveragingMethod='exponential',ForgettingFactor=0.99, ...
%             YLimits=[-30 10],ShowLegend=true);
% 
%     spectrumScope([cos_table.', COSTAB.']);
%%
next_phase_correction = 0;
pi_int = 1608;

for i = 1:SIZE
    dds_cos1(i) = i2(index);
    dds_sin1(i) = q2(index);

%     if (next_phase_correction < -pi_int)
%         next_phase_correction = next_phase_correction + 2*pi_int;
%     end
% 
%     next_phase_correction = next_phase_correction - step;
%     index = abs(next_phase_correction);

    index = index+step;

    if index > SIZE
        index = index + step - SIZE;
    end
end

radians_per_sample = 0.0314;
time_base=0:SIZE-1;
correction_signal=exp(-j*(radians_per_sample)*time_base);

nn = (1:SIZE).';
figure(12);
a1 = dds_cos1;
a1 = a1;
% a11 = time_base(1:end-1).';
% a111 = a1.*a11;
% a2 = real(correction_signal*2^11);
a2 = real(correction_signal);
a2 = a2.';
as = angle(correction_signal);
ar = angle(complex(dds_cos1,dds_sin1));
% plot(nn, ar, nn, as);
plot(nn, a1, nn, a2);

% figure(14)
% n11 = (1:256).';
% cosi = importdata('test_signals\rotate_mult_i.txt');
% plot(n11, cosi/2048, n11, dds_cos1(1:256))
