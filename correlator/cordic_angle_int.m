% CORDIC_For_Dummies.pdf
% Doc 51 - Atan2 Cordic Algorithm.pdf
% https://dspguru.com/dsp/faqs/cordic/
% https://openofdm.readthedocs.io/en/latest/verilog.html

%%
function angle_int = cordic_angle_int(i, q, iterations)
   
% i = int32(-81.2793*2^11);
% q = int32(-4.089*2^11);
% iter = 8;
% cordic_matlab = cordicangle(complex(i,q),iter)*180/pi;

    % create atan table in int
    % atan_table
    % atan(2^0) = 0.7854 = 0.7854 * 180/pi = 45 degrees * 256 = 11520
    % atan(2^-1) = 0.4636 = 0.4636 * 180/pi = 26.5651 degrees * 256 = 6801
    % atan(2^-2) = 0.2450 = 0.2450 * 180/pi = 14.0362 degrees * 256 = 3594
    % atan(2^-3) = 0.1244 = 0.1244 * 180/pi = 7.125 degrees * 256 = 1825
    % atan(2^-4) = 0.0624 = 0.0624 * 180/pi = 3.5763 degrees * 256 = 916
    % atan(2^-5) = 0.0312 = 0.0312 * 180/pi = 1.7899 degrees * 256 = 459

for j = 0:14
    atan_table_int(j+1) = int32(ceil((atan(2^-j)*180/pi)*256)); 
end
%% https://dspguru.com/dsp/faqs/cordic/   3.2.2

  i_tmp = i;

  if (i < 0)
        if (q > 0)
            % поворот на -90 градусов
            i = q;         % I = Q
            q = -i_tmp;    % Q = - I
            angle_int = int32(90*256); % начальный вектор имеет фазу в 90 градусов
        else
            % поворот на 90 градусов 
            i = -q;        % I = - Q;
            q = i_tmp;     % Q = I
            angle_int = int32(-90*256); % начальный вектор имеет фазу в -90 градусов
        end
  else
        angle_int = int32(0);
  end
%% https://dspguru.com/dsp/faqs/cordic/   3.2.3

  for k = 0:iterations-1
  % против часовой стрелки
    if (q >= 0)
        in = i + (bitshift(q,-k,'int32'));
        qn = q - (bitshift(i,-k,'int32'));
        angle_int = angle_int + atan_table_int(k+1);
    else
        % по часовой стрелке
        in = i - (bitshift(q,-k,'int32'));
        qn = q + (bitshift(i,-k,'int32'));
        angle_int = angle_int - atan_table_int(k+1);
    end

    i = in;
    q = qn;

  end