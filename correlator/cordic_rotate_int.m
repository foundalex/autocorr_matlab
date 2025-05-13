
clear all;
phase_rad = 0.0503;
phase_deg = floor((phase_rad*180/pi)*256);

niters = 5;

i = int32(0.0469 * 2^11);
q = int32(0.0877 * 2^11);
signal = 0.0469+0.0877i;

a = cordicrotate(phase_rad, signal, niters);

for j = 0:14
    atan_table_int(j+1) = int32(ceil((atan(2^-j)*180/pi)*256)); 
end

AnGain = floor(256 * prod(sqrt(1+2.^(-2*(0:(niters-1))))));

   if (phase_deg > 90*256)
        xt = 0;
        yt = 256;
        angle1 = 90*256;
   elseif (phase_deg < -90*256)
        xt = 0;
        yt = -256;
        angle1 = -90*256;
   else
        xt = 256;
        yt = 0;
        angle1 = 0;
   end

 for k = 0:niters-1
        % против часовой стрелки
        if (phase_deg - angle1 < 0)
            xn_t = xt + (bitshift(yt,-k,'int32'));
            yn_t = yt - (bitshift(xt,-k,'int32'));
            angle1 = angle1 - atan_table_int(k+1);
        % по часовой стрелке    
        else
            xn_t = xt - (bitshift(yt,-k,'int32'));
            yn_t = yt + (bitshift(xt,-k,'int32'));
            angle1 = angle1 + atan_table_int(k+1);
        end

        xt = xn_t;
        yt = yn_t;

        p_cos = xt / AnGain;
        p_sin = yt / AnGain;

        a1 = signal * complex(p_cos,p_sin);
    end
