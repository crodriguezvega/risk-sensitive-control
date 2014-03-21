clear;

test = riskSensitiveHJB2factor();   
test.t = linspace(1 - 1/100, 1, 2);    
test.x1 = linspace(-10, 10, 50);
test.x2 = linspace(-5, 5, 50);  
test.v = 1000;
test.numberOfControls = 50;
test.theta = 2;
test.b = [0.3
          0.7];
test.B = [0.5 -0.3
          -0.3 0.2];
test.lambda = [0.2 -0.2
               -0.2 0.5];
test.a0 = 0.05;
test.A0 = [-0.3
           -0.4];
test.a = [0.2
          0.4];
test.A = [0.4 -0.3
          -0.2 0.5];
test.sigma = [0.3 -0.2
              -0.2 0.4];
test.z1 = linspace(-5, 5, 20);
test.z2 = linspace(-5, 5, 20);
test.gamma_min = -0.9;
test.gamma_max = 10;
test.gamma = (test.gamma_max - test.gamma_min) ...
             .*rand(length(test.z1), length(test.z2), 2) + test.gamma_min;
test.R = [1
          2];
test.intensity = 5;
test.mean = [0.5    
             0.8]; 
test.covariance = [0.8 0.5
                   0.5 0.8];

% Find solution to HJB problem
[V, h] = test.solve();

axis tight
[X, Y] = meshgrid(test.x2, test.x1); 
surf(X, Y, V(:, :, 1));
xlabel('Factor nr. 2');
ylabel('Factor nr. 1');
title('Value function');

figure;
surf(X, Y, h(:, :, 1, 1));
xlabel('Factor nr. 2');
ylabel('Factor nr. 1');
title('Optimal control for asset nr. 1');

figure;
surf(X, Y, h(:, :, 1, 2));
xlabel('Factor nr. 2');
ylabel('Factor nr. 1');
title('Optimal control for asset nr. 2');

figure;
surf(X, Y, ones(length(X), length(Y)) - (h(:, :, 1, 1) + h(:, :, 1, 2)));
xlabel('Factor nr. 2');
ylabel('Factor nr. 1');
title('Weight for money market account');