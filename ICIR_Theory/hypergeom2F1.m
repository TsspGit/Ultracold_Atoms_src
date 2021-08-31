function res=hypergeom2F1(a,b,c,z)
% Faster hyper geometric function for real arguments.
tol=eps;
if z<eps || max([imag(a),imag(b),imag(c),imag(z)])>eps || mod(a-b,1)==0 || mod(c-a-b,1)==0
    disp('Parameter beyond the scope of this function');
    res=NaN;
end

if (z<2) && (z>1)
    fac2=-2*b;
    if b>-100
        fac1=sqrt(pi)*gamma(1+b)/gamma(0.5+b);
    else
        fac1=-sqrt(pi)*(0.5+b)*sin(pi*(0.5+b))/sin(pi*b);
        z=-0.5-b;
        gamma1=(z-0.5)*log(z)-z+0.5*log(2*pi)+1/(12*z)-1/(360*z^3)+1/(1260*z^5);
        z=1-b;
        gamma2=(z-0.5)*log(z)-z+0.5*log(2*pi)+1/(12*z)-1/(360*z^3)+1/(1260*z^5);
        fac1=fac1*exp(gamma1-gamma2);
    end
    w=1-1/z;
    res1=hypergeom2F1_Horchler(a,a-c+1,a+b-c+1,w,tol)*(1-w)^a*fac1;
    res2=hypergeom2F1_Horchler(1-b,c-b,c-b-a+1,w,tol)*(1-w)^a*w^(c-a-b)*exp(1i*pi*(c-a-b))*fac2;
    res=res1+res2;
else
    fac1=2*b/(2*b-1);
    if b>-100
        fac2=gamma(1+b)*gamma(0.5-b)/sqrt(pi);
    else
        z=0.5-b;
        gamma1=(z-0.5)*log(z)-z+0.5*log(2*pi)+1/(12*z)-1/(360*z^3)+1/(1260*z^5);
        z=1-b;
        gamma2=(z-0.5)*log(z)-z+0.5*log(2*pi)+1/(12*z)-1/(360*z^3)+1/(1260*z^5);
        fac2=b*sqrt(pi)*exp(gamma1-gamma2)/sin(b*pi);
    end
    if (z<0) && (z>=-1)
        w=z/(z-1);
    elseif (z<=0.5) && (z>=0)
        w=z;
    elseif (z<=1) && (z>0.5)
        w=1-z;
    elseif z>2
        w=1/z
    end
    %disp(w)
    res1=hypergeom2F1_Horchler(a,a-c+1,a-b+1,w,tol)*abs(w)^a*fac1*exp(1i*pi*a);
    res2=hypergeom2F1_Horchler(b-c+1,b,b-a+1,w,tol)*abs(w)^b*fac2*exp(1i*pi*b);
    res=(res1+res2)';
end

if isnan(res)||isinf(res)
    disp(['Hypergeom2F1 goes wrong when [a,b,c,z]=',num2str([a,b,c,z])]);
    res=NaN;
end

end

function h=hypergeom2F1_Horchler(a,b,c,z,tol)
%HYPERGEOM2F1  Gaussian or ordinary hypergeometric function for ABS(Z) < 1
%   H = HYPERGEOM2F1(A,B,C,Z) returns the hypergeometric function 2F1 for scalar
%   parameters A, B, C, and complex inputs Z. A finite number of terms from the
%   infinite summation representation of the function are used until the
%   absolute tolerance EPS is met.
%
%   H = HYPERGEOM2F1(...,TOL) specifies the absolute tolerance used to terminate
%   the summation of terms.
%   
%   Note:
%       Unless, C = A or C = B, if C is an integer <= 0, NaN is returned.
%       Additionally, the simple method of computation used can be very
%       inaccurate when ABS(Z) is close to 1 for some parameter combinations.
%
%   See also HYPERGEOM.

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-12-13
%   Revision: 1.0, 5-19-14


% Check four required inputs
if isempty(a) || ~isscalar(a) || ~isfloat(a) || ~isreal(a) || ~isfinite(a)
    error('hypergeom2F1:AInvalid',...
          'A must be a non-empty finite real floating-point scalar.');
end
if isempty(b) || ~isscalar(b) || ~isfloat(b) || ~isreal(b) || ~isfinite(b)
    error('hypergeom2F1:BInvalid',...
          'B must be a non-empty finite real floating-point scalar.');
end
if isempty(c) || ~isscalar(c) || ~isfloat(c) || ~isreal(c) || ~isfinite(c)
    error('hypergeom2F1:CInvalid',...
          'C must be a non-empty finite real floating-point scalar.');
end
if ~isfloat(z) || ~all(isfinite(z))
    error('hypergeom2F1:ZInvalid','Z must be a finite floating-point array.');
end
if any(abs(z)>=1)
    error('hypergeom2F1:ZUndefined',...
         ['The standard Gaussian hypergeometric function is only defined ' ...
          'for ABS(Z) < 1']);
end

% Calculate Gaussian hypergeometric function via summation of terms
if isempty(z)
	h = z;
else
    % Set relative tolerance and maximum number of iterations
    dtype = superiorfloat(a,b,c,z);
    if nargin < 5
        tol = eps(dtype);
    end
    itermax = 2^15;
    
    h = zeros(size(z),dtype);
    for j = 1:numel(z)
        Z = z(j);
        
        if a == 1
            if b == c
                yi = Z;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*Z;
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            else
                yi = b*Z/c;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*(b+i)*Z/(c+i);
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            end
        elseif b == 1
            if a == c
                yi = Z;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*Z;
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            else
                yi = a*Z/c;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*(a+i)*Z/(c+i);
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            end
        elseif a == c
            yi = b*Z;
            y = 1+yi;
            for i = 1:itermax
                yi = yi*(b+i)*Z/(i+1);
                y = y+yi;
                if abs(yi) < tol
                    break;
                end
            end
        elseif b == c
            yi = a*Z;
            y = 1+yi;
            for i = 1:itermax
                yi = yi*(a+i)*Z/(i+1);
                y = y+yi;
                if abs(yi) < tol
                    break;
                end
            end
        else 
            yi = a*b*Z/c;
            y = 1+yi;
            for i = 1:itermax
                yi = yi*(a+i)*(b+i)*Z/((i+1)*(c+i));
                y = y+yi;
                if abs(yi) < tol
                    break;
                end
            end
        end
        
        h(j) = y;
    end
end
end