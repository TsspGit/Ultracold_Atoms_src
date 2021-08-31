function [aq]=adquad2f1(a,b,c,z,h,tol1,tol2)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) for             %
% parameters a,b and variable z by applying the method of adaptive        %
% quadrature as detailed in Appendix H.1.                                 %

%-------------------------------------------------------------------------%
% Input:  a=Parameter value a                                             %
%         b=Parameter value b                                             %
%         c=Parameter value c                                             %
%         z=Variable value z                                              %
%         h=Initial step-size                                             %
%         tol1=Tolerance as defined in Appendix G.5                       %
%         tol2=Tolerance as defined in Appendix G.5                       %
% Output: aq=Computed value of 2F1(a,b;c;z)                               %
%-------------------------------------------------------------------------%

% Define integrand
f=@(y) (1-z*y)^(-a)*y^(b-1)*(1-y)^(c-b-1);

% Start at left of integral (4.8)
x=0;
accumulator=0;
% If reach right of interval [0,1], don't go past it
while x<1
    if x+h>1
        h=1-x;
    end
    % Two quadrature methods used, composite trapezoidal and composite
    % Simpson's
    quad1=h/2*(f(x)+f(x+h));
    quad2=h/6*(f(x)+4*f(x+h/2)+f(x+h));
    % Decrease step-size if difference between two numerical approximations
    % is too great
    if abs(quad1-quad2)>tol1
        h=h/2;
    else
        accumulator=accumulator+quad2;
        x=x+h;
        % Increase strp-size if difference is small
        if abs(quad1-quad2)<tol2
            h=2*h;
        end
    end
end

% Return solution to (4.8)
aq=accumulator*gamma(c)/gamma(b)/gamma(c-b);