
a = 100;
b = -100;
n = 100;
niters = 1;


phase = 81.2793 + 4.089i;
cordic_user = cordic_angle(phase,niters)*180/pi;
cordic_matlab = cordicangle(phase,niters)*180/pi;

for k = 1:n
    i = a + (b-a).*rand(1,1);
    q = a + (b-a).*rand(1,1);
    iter = randi([1 8],1,1);
    phase = complex(i,q);
    cordic_user = cordic_angle(phase,iter);
    cordic_matlab = cordicangle(phase,iter);
    error_cordic = cordic_matlab - cordic_user;
    if (error_cordic > 1^-14)
        X = [num2str(k), "Test failed"];
        disp(X);
    end
end