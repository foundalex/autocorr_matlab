function [dds_cos, dds_sin] = dds_deg_int(angle_int,n)
% 1/4 period of sine
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
% fo = 100000;
% fs = 20000000;
% ang_rad = (-fo*2*pi*16)/fs;
% ang_rad_w = ang_rad/16;
% n = 1024;
% ang_rad = -7376;
% ang_rad_w = ang_rad/16
ang_rad = angle_int*pi/180/256*16;
ang_rad_w = floor((angle_int*2)*pi/180);
pi_int = 1608;
%%
index = 0;
next_phase_correction = 0;

for i = 1:n

    if (ang_rad < 0)
        next_phase_correction = next_phase_correction + round(ang_rad_w*512);
        if (next_phase_correction < -pi_int)
            next_phase_correction = next_phase_correction + 2*pi_int;
        end
    end

    if abs(next_phase_correction) <= pi_int/4
        if (next_phase_correction <= 0)
            quadrant = 4;
        else
            quadrant = 0;
        end
        index = abs(next_phase_correction);
        if (index == 0)
            index = abs(round(ang_rad_w*512));
        end
    elseif abs(next_phase_correction) <= pi_int/2
        if (next_phase_correction < 0)
            quadrant = 5;
        else
            quadrant = 1;
        end
        index = pi_int/2 - abs(next_phase_correction);
    elseif abs(next_phase_correction) <= pi_int*3/4
        if (next_phase_correction < 0)
            quadrant = 6;
        else
            quadrant = 2;
        end
          index = abs(next_phase_correction) - pi_int/2;
    else
        if (next_phase_correction < 0)
            quadrant = 7;
        else
            quadrant = 3;
        end
        index = pi_int - abs(next_phase_correction);
    end
    
    sin1(i) = str_i_array(index+1);
    cos1(i) = str_q_array(index+1);

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
end

dds_cos = int16(rot_i).';
dds_sin = int16(rot_q).';