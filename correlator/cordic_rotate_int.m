
function [p_cos, p_sin] = cordic_rotate_int(phase_deg, gain, niters)

scale = 2^28;
pi_int = pi*scale;

    for j = 0:niters
        atan_table_int(j+1) = int32(ceil(atan(2^-j)*scale)); 
    end

    if (phase_deg > pi_int/2)
        xt = 0;
        yt = scale;
        angle1 = pi_int/2;
    elseif (phase_deg < -pi_int/2)
        xt = 0;
        yt = -scale;
        angle1 = -pi_int/2;
    else
        xt = scale;
        yt = 0;
        angle1 = 0;
    end

    for k = 0:niters-1
        % против часовой стрелки
        if (phase_deg - angle1 < 0)
            xn_t = xt + (bitshift(yt,-k,'int64'));
            yn_t = yt - (bitshift(xt,-k,'int64'));
            angle1 = angle1 - atan_table_int(k+1);
        % по часовой стрелке    
        else
            xn_t = xt - (bitshift(yt,-k,'int64'));
            yn_t = yt + (bitshift(xt,-k,'int64'));
            angle1 = angle1 + atan_table_int(k+1);
        end

        xt = xn_t;
        yt = yn_t;

        p_cos = xt / gain;
        p_sin = yt / gain;
 
    end
