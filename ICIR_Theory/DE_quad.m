%% author: Tomás Sánchez Sánchez-Pastor 31/08/21
clear all
clc
%% Function to integrate
x  = linspace(0, 1000, 100000);
h  = x(end) - x(end-1); % step
%% Doubly exponential quadrature || tanh-sinh quad
% v1
I  = 0;
for k=-x(end):x(end)
   xk = exp(pi/2 * sinh(k*h));
   wk = pi/2*h*cosh(k*h)*xk;
   I = I + get_f(xk) * wk;
end
disp(I)

% v2
kmax = 200000;
h = 2e-5;
I = 0;
for k=-kmax:kmax
    xk = exp(pi/2 * sinh(k*h));
    wk = pi/2*h*cosh(k*h)*xk;
    I = I + get_f(xk) * wk;
end
disp(I)

function f = get_f(x)
    f = exp(-x^2);
end