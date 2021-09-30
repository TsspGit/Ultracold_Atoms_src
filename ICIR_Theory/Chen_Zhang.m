%% autor: Tomás Sánchez Sánchez-Pastor, 29-07-21
clc
clear all
%% Parameters:
etax  = 1;
etay  = 1;
etaz  = 1;
nx     = 100;
ny     = 100;
Eo     = 1/2*(etax + etay + etaz);
E      = linspace(-7.5, 12, 500);

%% Integration of A3D
a3D = zeros(1, length(E));
B1_3D = zeros(1, length(E));
B2_3D = zeros(1, length(E));
A3D_int = zeros(1, length(E));
for k = 1:length(E)
    if E(k) < Eo
        Lambda = 100;
        a = 1e-6;
        [beta, w] = lgwt(500, a, Lambda);
        I3D_int(k) = sum(I3D_f(etax, etay, nx, ny, E(k), beta).*w);
    elseif E(k) >= Eo
        Lambda = 0.1/(E(k) - Eo);
        a = 1e-6;
        [beta, w] = lgwt(500, a, Lambda);
        I3D_int(k) = sum(I3D_f(etax, etay, nx, ny, E(k), beta).*w);
    end
    B1_3D(k) = B1_3D_f(nx, ny, etax, etay, E(k), Lambda);
    B2_3D(k) = B2_3D_f(nx, ny, etax, etay, E(k), Lambda);
    J3D = sqrt(2)*4*pi.*(W3D_f(etax, etay, nx, ny, E(k)) + I3D_int(k));
    a3D(k) = 1./J3D;
end

%% Plot
% B1_3D
figure
grid on
plot(1./B1_3D, E, 'b-', 'linewidth', 2)
xlabel('$d_y/B1_{3D}$', 'Interpreter', 'latex')
ylabel('$E/(\hbar \omega_z)$', 'Interpreter', 'latex')
xlim([-10, 10])

% B2_3D
figure
grid on
plot(1./B2_3D, E, 'b-', 'linewidth', 2)
xlabel('$d_y/B2_{3D}$', 'Interpreter', 'latex')
ylabel('$E/(\hbar \omega_z)$', 'Interpreter', 'latex')
xlim([-10, 10])
% E(a3D)
figure
grid on
plot(a3D, E, 'b-', 'linewidth', 2)
xlabel('$a_{3D}/d_y$', 'Interpreter', 'latex')
ylabel('$E/(\hbar \omega_z)$', 'Interpreter', 'latex')
xlim([-10, 10])
ylim([-7.5, 12])
saveas(gcf, 'as_E_isotropic.png')
%exit;

%% Functions
function en = en_f(n, eta)
% n:   level
% eta: anisotropy
en = eta.*(n + 1./2);
end

function A3D = A3D_f(etax, etay, E, beta)
        A3D = 1/(4*pi)^(3/2) * (-exp(beta.*E).* sqrt(etax*etay./(sinh(beta.*etax).*sinh(beta.*etay)...
        .*sinh(beta))) + beta.^(-3/2));
end

function I3D = I3D_f(etax, etay, nx, ny, E, beta)
        suma = 0;
        for i=0:2:nx
            for j=0:2:ny
                if E >= en_f(i, etax) + en_f(j, etay) + 1/2
                    suma = suma + 2^(i+j-1/2) * exp(beta*(E-en_f(i, etax) - en_f(j, etay)))./(gamma(1/2-i/2).^2 .* gamma(1/2-j/2).^2 ...
                     .*gamma(i+1).*gamma(j+1));
                end
            end
        end
        I3D = A3D_f(etax, etay, E, beta) + sqrt(pi*etax*etay./(8.*sinh(beta))).*suma;
end

function W3D = W3D_f(etax, etay, nx, ny, E)
    suma = 0;
    for i=0:nx
        for j=0:ny
            if ((mod(i, 2) == 0) && (mod(j, 2) == 0) && (E >= en_f(i, etax) + en_f(j, etay) + 1/2))
                suma = suma + 2^(i+j-1) * gamma(1/4 - (E - en_f(i, etax) - en_f(j, etay))/2)./(gamma(1/2-i/2).^2 .* gamma(1/2-j/2).^2 ...
                 .*gamma(i+1).*gamma(j+1).*gamma(3/4 - (E - en_f(i, etax) - en_f(j, etay))/2));
            end
        end
    end
    W3D = -pi/2 * sqrt(etax*etay/2) * suma;
end

function B1_3D = B1_3D_f(nx, ny, etax, etay, E, Lambda)
% nj:     Levels
% etaj:   Anisotropy wj./wy
% E:      Energy
% Eo:     Free energy
% Lambda: Integration limit
sum = 0;
for i = 0:2:nx
    for j = 0:2:ny
        if E >= en_f(i, etax) + en_f(j, etay) + 1/2
             sum = sum + (2^(i + j - 5/2).*sqrt(pi*etax*etay) * gamma(1/4 - (E - en_f(i, etax) - en_f(j, etay))/2) ...
                 .*exp((E - en_f(i, etax) - en_f(j, etay) - 3/2).*Lambda))./( gamma(1/2-i/2).^2 .* gamma(1/2-j/2).^2 ...
                 .*gamma(i+1).*gamma(j+1).*gamma(5/4 - (E - en_f(i, etax) - en_f(j, etay))/2)) .* sqrt(exp(2.*Lambda)...
                 - 1) .* hypergeom2F1(1, 3./4 - (E - en_f(i, etax) - en_f(j, etay))./2, 5./4 - (E - en_f(i, etax) - ...
                 en_f(j, etay))./2, exp(-2.*Lambda));

        end
    end
end
B1_3D = (-1) .* sum;
end

function B2_3D = B2_3D_f(nx, ny, etax, etay, E, Lambda)
% nj:     Levels
% etaj:   Anisotropy wj./wy
% E:      Energy
% Eo:     Free energy
% Lambda: Integration limit
sum = 0;
for i = 0:1:nx
    for j = 0:1:ny
        if ((mod(i, 2)==1) && (mod(j, 2) == 1) || ((mod(i, 2)==0) && (mod(j, 2) == 0) && (E<en_f(i, etax) + en_f(j, etay) + 1/2)))
            sum = sum + ( 2^(i+j-1).*exp((E - en_f(i, etax) - en_f(j, etay) - 2).*Lambda).*(exp(2.*Lambda)-1)...
                .* hypergeom2F1(1, 3/4 - (E - en_f(i, etax) - en_f(j, etay))/2, 5/4 - (E - en_f(i, etax) - ...
                en_f(j, etay))/2, exp(-2.*Lambda))) ./ ( gamma(1/2-i/2).^2 .* gamma(1/2-j/2).^2 ...
                .*gamma(i+1).*gamma(j+1).*(2.*(E-en_f(i, etax)-en_f(j, etay))-1) );
        end
    end
end
B2_3D = sqrt(pi.*etax.*etay./sinh(Lambda)) .* sum;
end

function [x,w]=lgwt(N,a,b)
% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;
xu=linspace(-1,1,N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
% Derivative of LGVM
Lp=zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end
% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      
% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end