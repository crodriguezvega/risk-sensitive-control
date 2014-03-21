clear;

test = riskSensitiveHJB1factor();   
test.t = linspace(0, 1, 100);    
test.x1 = linspace(-10, 10, 50);
test.v = 1000;
test.numberOfControls = 500;
test.theta = 2;
test.b = 0.7;
test.B = 0.5;
test.lambda = 0.2;
test.a0 = 0.05;
test.A0 = 0.02;
test.a = 0.4;
test.A = 0.8;
test.sigma = 0.3;
test.z = linspace(-5, 5, 100);
test.gamma_min = -0.9;
test.gamma_max = 10;
test.gamma = (test.gamma_max - test.gamma_min) ...
             .*rand(length(test.z), 1) + test.gamma_min;
test.R = 2;
test.intensity = 5;
test.mean = 1; 
test.variance = 0.8;

% Find solution to HJB problem
[V, h] = test.solve();

axis tight
[t, X] = meshgrid(test.t, test.x1); 
surf(t, X, V);
xlabel('Time');
ylabel('Factor');
title('Value function');

figure
[t, X] = meshgrid(test.t(1:end -1), test.x1); 
surf(t, X, h);
xlabel('Time');
ylabel('Factor');
title('Optimal control');