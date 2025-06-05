function [dds_cos, dds_sin] = dds_int(angle_rad_int, n)
% 1/4 period of sine import from verilog
lines = readlines("test_signals\rot_lut.txt");
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
%%
% SCALE = 2048;
% MAX3 = round(pi/4*SCALE);
% SIZE1 = 2^round(log2(MAX3));
% for k = 1:SIZE1
%     key2(k) = k/MAX3*pi/4;
%     i_I(k) = (round(cos(key2(k))*2^11));
%     i_Q(k) = (round(-sin(key2(k))*2^11));
% end
% 
% pi_int = (round(pi*2^11));
% pi_int_2 = (2*(round(pi*2^11)));
% 
% pi_int_div_4 = (round((pi*2^11)/4));
% pi_int_div_2 = (round((pi*2^11)/2));
% pi_int_div_3_4 = (round((pi*2^11)*3/4));

% r1 = str_i_array(1:512)-double(i_I2(1:512));

%%
% next_phase_correction = 0;
% for i = 2:n
%         next_phase_correction = next_phase_correction + angle_rad_int;
%         if (next_phase_correction < -SIZE1)
%             next_phase_correction = next_phase_correction + SIZE1;
%         elseif (next_phase_correction > SIZE1)
%             next_phase_correction = next_phase_correction - SIZE1;
%         end
%     index = (abs(next_phase_correction));
%     dds_cos(i) = i_I(index);
%     dds_sin(i) = i_Q(index);
% end
% 
% dds_cos = dds_cos.';
% dds_sin = dds_sin.';

%% verilog model
pi_int = int16(1608);
pi_int_2 = int16(3216);

pi_int_div_4 = 402;
pi_int_div_2 = 804;
pi_int_div_3_4 = 1206;

index = 0;
next_phase_correction = int32(0);
dds_cos(1) = str_i_array(index+1);
dds_sin(1) = str_q_array(index+1);
for i = 2:n

    if (angle_rad_int < 0)
        next_phase_correction = next_phase_correction + int32(angle_rad_int);
        if (next_phase_correction < int32(-pi_int))
            next_phase_correction = next_phase_correction + int32(pi_int_2);
        elseif (next_phase_correction > int32(pi_int))
            next_phase_correction = next_phase_correction - int32(pi_int_2);
        end
    end

    if abs(next_phase_correction) <= pi_int_div_4 %(pi_int/4)
        if (next_phase_correction <= 0)
            quadrant = 4;
        else
            quadrant = 0;
        end
        index = (abs(next_phase_correction));
    elseif abs(next_phase_correction) <= pi_int_div_2 %pi_int/2
        if (next_phase_correction < 0)
            quadrant = 5;
        else
            quadrant = 1;
        end
        index = (int32(pi_int_div_2) - abs(next_phase_correction));
    elseif abs(next_phase_correction) <= pi_int_div_3_4 %pi_int*3/4
        if (next_phase_correction < 0)
            quadrant = 6;
        else
            quadrant = 2;
        end
        index = (abs(next_phase_correction) - int32(pi_int_div_2));
    else
        if (next_phase_correction < 0)
            quadrant = 7;
        else
            quadrant = 3;
        end
        index = int32(pi_int) - abs(next_phase_correction);
    end
    
    sin1(i) = str_i_array((index)+1);
    cos1(i) = str_q_array((index)+1);

        if (quadrant == 0)
            rot_i(i) = sin1(i);
            rot_q(i) = cos1(i);
        elseif (quadrant == 1)
            rot_i(i) = cos1(i);
            rot_q(i) = sin1(i);
        elseif (quadrant == 2)
            rot_i(i) = -cos1(i);
            rot_q(i) = sin1(i);
        elseif (quadrant == 3)
            rot_i(i) = -sin1(i);
            rot_q(i) = cos1(i);
        elseif (quadrant == 4)
            rot_i(i) = sin1(i);
            rot_q(i) = -cos1(i);
        elseif (quadrant == 5)
            rot_i(i) = cos1(i);
            rot_q(i) = -sin1(i);
        elseif (quadrant == 6)
            rot_i(i) = -cos1(i);
            rot_q(i) = -sin1(i);
        elseif (quadrant == 7)
            rot_i(i) = -sin1(i);
            rot_q(i) = -cos1(i);
        end

%         quadrant_i(i) = quadrant*100;
end

dds_cos(2:n) = rot_i(2:n);
dds_sin(2:n) = rot_q(2:n);

