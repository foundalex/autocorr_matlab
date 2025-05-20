function [p_cos, p_sin] = cordic_rotate(phase_rad, gain, niters)


tan_table = [1; 0.5; 0.25; 0.125; 0.0625; 0.0312; 0.0156; 0.0078];

atan_table = [0.785398163397448; 0.463647609000806; 0.244978663126864; 0.124354994546761; 
              0.062418809995957; 0.031239833430268; 0.015623728620477; 0.007812341060101];

   if (phase_rad > pi/2)
        xt = 0;
        yt = 1;
        angle1 = pi/2;
   elseif (phase_rad < -pi/2)
        xt = 0;
        yt = -1;
        angle1 = -pi/2;
   else
        xt = 1;
        yt = 0;
        angle1 = 0;
   end


    for i = 1:niters
        % против часовой стрелки
        if (phase_rad - angle1 < 0)
            xn_t = xt + yt * tan_table(i);
            yn_t = yt - xt * tan_table(i);
            angle1 = angle1 - atan_table(i);
        % по часовой стрелке    
        else
            xn_t = xt - yt * tan_table(i);
            yn_t = yt + xt * tan_table(i);
            angle1 = angle1 + atan_table(i);
        end

        xt = xn_t;
        yt = yn_t;

        p_cos = xt / gain;
        p_sin = yt / gain;
% 
%         a1 = signal * complex(p_cos,p_sin);
    end

    