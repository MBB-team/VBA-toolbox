function   y = sigmoid(x)
% sigmoid mapping [see invsigmoid.m]
y=1./(1+exp(-x));
y(y<1e-4) = 1e-4;
y(y>(1-1e-4)) = 1-1e-4;
