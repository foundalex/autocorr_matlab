% CORDIC_For_Dummies.pdf
% Doc 51 - Atan2 Cordic Algorithm.pdf
% https://dspguru.com/dsp/faqs/cordic/

function angle1 = cordic_angle(phase, iterations)

   tan_table = [1; 0.5; 0.25; 0.125; 0.0625; 0.0312; 0.0156; 0.0078];

    % atan_table
    % atan(2^0) = 0.7854 = 0.7854 * 180/pi = 45 degrees
    % atan(2^-1) = 0.4636 = 0.4636 * 180/pi = 26.5651 degrees
    % atan(2^-2) = 0.2450 = 0.2450 * 180/pi = 14.0362 degrees
    % atan(2^-3) = 0.1244 = 0.1244 * 180/pi = 7.125 degrees
    % atan(2^-4) = 0.0624 = 0.0624 * 180/pi = 3.5763 degrees
    % atan(2^-5) = 0.0312 = 0.0312 * 180/pi = 1.7899 degrees
    
   atan_table = [0.7854; 0.4636; 0.2450; 0.1244; 
                 0.0624; 0.0312; 0.0156; 0.0078];

   xt = real(phase);
   yt = imag(phase);

   x_tmp = xt;
       
  % https://dspguru.com/dsp/faqs/cordic/   3.2.2
  
  if (xt < 0)
        if (yt > 0)
            % поворот на -90 градусов
            xt = yt;        % I = Q
            yt = -x_tmp;    % Q = - I
            angle1 = pi/2; % начальный вектор имеет фазу в 90 градусов
        else
            % поворот на 90 градусов 
            xt = -yt;       % I = - Q;
            yt = x_tmp;     % Q = I
            angle1 = -pi/2; % начальный вектор имеет фазу в -90 градусов
        end
  else
        angle1 = 0;
  end
%%   
    % https://dspguru.com/dsp/faqs/cordic/   3.2.3

    for i = 1:iterations
    % против часовой стрелки
        if (yt >= 0)
            xn_t = xt + yt * tan_table(i); 
            yn_t = yt - xt * tan_table(i);
            angle1 = angle1 + atan_table(i);
        else
            % по часовой стрелке
            xn_t = xt - yt * tan_table(i);
            yn_t = yt + xt * tan_table(i);
            angle1 = angle1 - atan_table(i);
        end

        xt = xn_t;
        yt = yn_t;

    end