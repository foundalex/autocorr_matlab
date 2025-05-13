
clear all;

a = int32(2^11);
b = int32(-2^11);
n = 100;

% 
% phase = -81.2793 + -4.089i;
% cordic_user = cordic_angle(phase,8)*180/pi;
% cordic_matlab = cordicangle(phase,8)*180/pi;

for k = 1:n
    i = a + (b-a).*rand(1,1);
    q = a + (b-a).*rand(1,1);
    iter = randi([1 8],1,1);
    phase = complex(i,q);
    cordic_user = cordic_angle_int(i, q, iter)/256;
    cordic_matlab = cordicangle(phase,iter)*180/pi;
    error_cordic = cordic_matlab - cordic_user;
    if (abs(error_cordic) > 0.5)
        X = [num2str(k), "Test failed"];
        disp(X);
    end
end