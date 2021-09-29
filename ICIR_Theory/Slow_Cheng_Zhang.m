%% autor: Tomás Sánchez Sánchez-Pastor, 15-09-21
clc
clear all
%% Parameters:
etax  = 1;
etay  = 1;
etaz  = 1;
nx     = 100;
ny     = 100;
Eo     = 1/2*(etax + etay + etaz);
E      = linspace(-7.5, 1.4, 600);
beta   = linspace(1e-5, 10, 1e3);
%% J3D
a3D = zeros(1, length(E));
for k=1:length(E)
    int_I3D = trapz(I3D_f(etax, etay, nx, ny, E(k), beta));
    J3D = 4*pi*(W3D_f(etax, etay, nx, ny, E(k)) + int_I3D);
    a3D(k) = 1/J3D;
end
%% Plot
figure
grid on
plot(a3D, E, 'b-', 'linewidth', 2)
xlabel('$a_{3D}/d_y$', 'Interpreter', 'latex')
ylabel('$E/(\hbar \omega_z)$', 'Interpreter', 'latex')
%xlim([-10, 10])
ylim([-7.5, 15])


%% Functions
function en = en_f(n, eta)
% n:   level
% eta: anisotropy
en = eta.*(n + 1./2);
end

function W3D = W3D_f(eta_x, eta_y, nx, ny, E)
% nj:     Levels
% etaj:   Anisotropy wj./wy
% E:      Energy
sum = 0;
for i=0:2:nx
    for j=0:2:ny
        if (en_f(i, eta_x) + en_f(j, eta_y) + 1/2) <= E
        sum = sum + 2^(i+j-1) * gamma(1/4 - (E-en_f(i, eta_x) - en_f(j, eta_y))/2) / (gamma(1/2-i/2)^2 ...
            * gamma(1/2-j/2)^2 * gamma(i+1) * gamma(j+1) * gamma(3/4 - (E-en_f(i, eta_x) - en_f(j, eta_y))/2));
        end
    end
end
W3D = -pi/2 * sqrt(eta_x*eta_y/2) * sum;
end

function I3D = I3D_f(eta_x, eta_y, nx, ny, E, beta)
% nj:     Levels
% etaj:   Anisotropy wj./wy
% E:      Energy
% beta:   X-axis
sum = 0;
for i=0:2:nx
    for j=0:2:ny
        sum = sum + 2^(i+j-1/2) * exp(beta*(E-en_f(i, eta_x) - en_f(j, eta_y))) ./ (gamma(1/2-i/2)^2 ...
            .* gamma(1/2-j/2)^2 .* gamma(i+1) .* gamma(j+1));
    end
end
I3D = -exp(beta.*E) .* sqrt(eta_x*eta_y./((4*pi)^3 .* sinh(eta_x*beta) .* sinh(eta_y*beta) .* sinh(beta))) + ...
    (1./(4*pi*beta)).^(3/2) .* sqrt(pi*eta_x*eta_y./(8*sinh(beta))) .* sum;
end