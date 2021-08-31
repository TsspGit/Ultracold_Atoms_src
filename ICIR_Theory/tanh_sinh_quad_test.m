%% autor: Tomás Sánchez Sánchez-Pastor 09/08/21
clear all
clc
%% Integral
N = 50000; 

%h = 2.0 / (N - 1)
h=0.0002; %(1/2^12)

% k ranges from -(N-1)/2 to +(N-1)/2
k = -1 * ((N - 1) / 2.0);
k_max  = ((N - 1) / 2.0);
thesum = 0;

% Loop across integration interval
actual_iter = 0;
while k < k_max + 1

    % Compute abscissa
    x_k = tanh(pi * 0.5 * sinh(k * h));

    % Compute weight
    numerator = 0.5 * h * pi * cosh(k * h);
    dcosh  = cosh(0.5 * pi * sinh(k * h));
    denominator = dcosh*dcosh;
    w_k =  numerator / denominator;

    thesum = thesum + w_k * func(x_k);
    myepsilon = abs(w_k * func(x_k));
    if mod(actual_iter, 10000) == 0 && actual_iter > k_max/2
        disp("Iteration =");
        disp(actual_iter);
        disp("epsilon =");
        disp(myepsilon);
    end
    k = k + 1;
    actual_iter = actual_iter + 1;
end
disp('Actual iterations = ');
disp(actual_iter);
disp("Integral = ");
disp(thesum)

function f = func(x)
    if x == 1 || x == -1
        f = 0;
    else
        f =  1 / sqrt(1 - x ^ 2);
    end
end