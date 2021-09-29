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
E      = linspace(-7.5, 15, 600);

%% Integration of A3D
a3D = zeros(1, length(E));
for k = 1:length(E)
    %disp(E(k));
    if E(k) < Eo
        Lambda = 4;
        beta   = linspace(1e-6, Lambda, 1e4);
        A3D    = 1/(4*pi)^(3/2) * (-exp(beta.*E(k)).* sqrt(etax*etay*etaz./(sinh(beta.*etax).*sinh(beta.*etay)...
            .*sinh(beta.*etaz))) + beta.^(-3/2));
        integral = trapz(A3D);
    elseif E(k) >= Eo
        Lambda = 8/(E(k) - Eo);
        beta   = linspace(1e-6, Lambda, 1e4);
        A3D    = 1/(4*pi)^(3/2) * (-exp(beta.*E(k)).* sqrt(etax*etay*etaz./(sinh(beta.*etax).*sinh(beta.*etay)...
            .*sinh(beta.*etaz))) + beta.^(-3/2));
        integral = trapz(A3D);
    end
    J3D = sqrt(2)*4*pi.*(integral + B1_3D_f(nx, ny, etax, etay, E(k), Lambda) + ...
        B2_3D_f(nx, ny, etax, etay, E(k), Lambda) + (1/(2*pi))^3/2 .* 1./sqrt(2*Lambda));
    a3D(k) = 1./J3D;
end

%% Plot
figure
grid on
plot(a3D, E, 'b-', 'linewidth', 2)
xlabel('$a_{3D}/d_y$', 'Interpreter', 'latex')
ylabel('$E/(\hbar \omega_z)$', 'Interpreter', 'latex')
%xlim([-10, 10])
ylim([-7.5, 15])
saveas(gcf, 'as_E_isotropic.png')
%exit;

%% Functions
function en = en_f(n, eta)
% n:   level
% eta: anisotropy
en = eta.*(n + 1./2);
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
        if mod(1 - (3/4 - (E - en_f(i, etax) - en_f(j, etay))/2), 1) ~= 0
            if ((mod(i, 2)==1) && (mod(j, 2) == 1) && (E>=en_f(i, etax) + en_f(j, etay) + 1/2))
                sum = sum + ( 2^(i+j-1).*exp((E - en_f(i, etax) - en_f(j, etay) - 2).*Lambda).*(exp(2.*Lambda)-1)...
                    .* hypergeom2F1(1, 3/4 - (E - en_f(i, etax) - en_f(j, etay))/2, 5/4 - (E - en_f(i, etax) - ...
                    en_f(j, etay))/2, exp(-2.*Lambda))) ./ ( gamma(1/2-i/2).^2 .* gamma(1/2-j/2).^2 ...
                    .*gamma(i+1).*gamma(j+1).*(2.*(E-en_f(i, etax)-en_f(j, etay))-1) );
            end
        end
    end
end
B2_3D = sqrt(pi.*etax.*etay./sinh(Lambda)) .* sum;
end