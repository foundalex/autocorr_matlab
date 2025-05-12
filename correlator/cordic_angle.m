  
 
function angle1 = cordic_angle(phase, iterations)
   %% cordic angle
   tan_table = [1; 0.5; 0.25; 0.125; 0.0625; 0.0312; 0.0156; 0.0078];
   atan_table = [0.7854; 0.4636; 0.2450; 0.1244; 0.0624; 0.0312; 0.0156; 0.0078];


%    phase = -81.2793 - 4.089i;
   tan_v = imag(phase)/real(phase);
   tan_v_degrees = tan_v * 180/pi;

%    interations = 1;

   xt = real(phase);
   yt = imag(phase);

   yt_tmp = yt; 
   cordicangle(phase,iterations);

   x_tmp = xt;
       
       if (xt < 0)
            if (yt > 0)
                xt = yt;
                yt = -x_tmp;
                angle1 = pi/2;
            else
                xt = -yt;
                yt = x_tmp;
                angle1 = -pi/2;
            end
       else
           angle1 = 0;
       end

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