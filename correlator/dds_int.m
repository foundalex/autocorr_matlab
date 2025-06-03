function [dds_cos, dds_sin] = dds_int(angle_rad_int, n)
% 1/4 period of sine import from verilog
% lines = readlines("test_signals\rot_lut.txt");
% phase_verilog = importdata("test_signals\phase_correct.txt");

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

% VAL = pi;
% SCALE = 2^11;
% ATAN_LUT_SCALE = 512;
% % MAX1 = round(pi/4*ATAN_LUT_SCALE);
% % SIZE = 2^ceil(log2(MAX1));

%% convert int sine to double
SCALE = 2^12;
MAX2 = round(2*pi*SCALE/8);
SIZE = 402*8;
for k = 1:SIZE
    key(k) = k/MAX2*2*pi;
    i_I(k) = int16(round(cos(key(k))*2^11));
    i_Q(k) = int16(floor(-sin(key(k))*2^11));
end

next_phase_correction = int32(0);
dds_cos(1) = i_I(1);
dds_sin(1) = 0;
%%
for i = 2:n
    if ((angle_rad_int) < 0)
        next_phase_correction = next_phase_correction - angle_rad_int;
        if (next_phase_correction < -SIZE)
            next_phase_correction = next_phase_correction + SIZE;
        elseif (next_phase_correction > SIZE)
            next_phase_correction = next_phase_correction - SIZE;
        end
    end
    index = (abs(next_phase_correction));
    dds_cos(i) = i_I(index);
    dds_sin(i) = i_Q(index);
end

dds_cos = dds_cos.';
dds_sin = dds_sin.';

r1 = str_i_array(1:512)-double(i_I(1:512));

%%
% pi_int = 1608;
% index = 0;
% next_phase_correction = 0;
% sin1(1) = str_i_array(index+1);
% cos1(1) = str_q_array(index+1);
% for i = 2:n
% 
%     if (angle_rad_int < 0)
%         next_phase_correction = next_phase_correction - (angle_rad_int);
%         if (next_phase_correction < -SIZE)
%             next_phase_correction = next_phase_correction + SIZE;
%         elseif (next_phase_correction > SIZE)
%             next_phase_correction = next_phase_correction - SIZE;
%         end
%     end
% 
%     if abs(next_phase_correction) <= SIZE/8
%         if (next_phase_correction <= 0)
%             quadrant = 4;
%         else
%             quadrant = 0;
%         end
%         index = floor(abs(next_phase_correction));
%     elseif abs(next_phase_correction) <= SIZE/2/2
%         if (next_phase_correction < 0)
%             quadrant = 5;
%         else
%             quadrant = 1;
%         end
%         index = floor(SIZE/2/2 - abs(next_phase_correction));
%     elseif abs(next_phase_correction) <= SIZE/2*3/4
%         if (next_phase_correction < 0)
%             quadrant = 6;
%         else
%             quadrant = 2;
%         end
%             index = floor(abs(next_phase_correction) - SIZE/2/2);
%     else
%         if (next_phase_correction < 0)
%             quadrant = 7;
%         else
%             quadrant = 3;
%         end
%         index = round(SIZE/2 - abs(next_phase_correction));
%     end
% 
%     phase_correction_matlab(i) = next_phase_correction;
%     
%     sin1(i) = str_i_array((index)+1);
%     cos1(i) = str_q_array((index)+1);
% 
%         if (quadrant == 0)
%             rot_i(i) = sin1(i);
%             rot_q(i) = cos1(i);
%         elseif (quadrant == 1)
%             rot_i(i) = cos1(i);
%             rot_q(i) = sin1(i);
%         elseif (quadrant == 2)
%             rot_i(i) = -cos1(i);
%             rot_q(i) = sin1(i);
%         elseif (quadrant == 3)
%             rot_i(i) = -sin1(i);
%             rot_q(i) = cos1(i);
%         elseif (quadrant == 4)
%             rot_i(i) = sin1(i);
%             rot_q(i) = -cos1(i);
%         elseif (quadrant == 5)
%             rot_i(i) = cos1(i);
%             rot_q(i) = -sin1(i);
%         elseif (quadrant == 6)
%             rot_i(i) = -cos1(i);
%             rot_q(i) = -sin1(i);
%         elseif (quadrant == 7)
%             rot_i(i) = -sin1(i);
%             rot_q(i) = -cos1(i);
%     end
% end
% 
% dds_cos1 = rot_i.';
% dds_sin1 = rot_q.';
% 
% error_phase = phase_verilog - phase_correction_matlab(1:256).';


