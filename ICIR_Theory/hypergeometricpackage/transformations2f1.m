function [h]=transformations2f1(ar,ai,br,bi,cr,ci,zr,zi,rho)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Applies transformation formulae detailed in Section 4.6.                %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         cr=Re(c)                                                        %
%         ci=Im(c)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         rho=Critical value for the absolute value of a transformed      %
%             variable at which point we allow the transformation formula %
%             to be used                                                  %
% Output: h=Computed value of 2F1(a,b;c;z)                                %
%-------------------------------------------------------------------------%

% Compute a,b,c,z in terms of ar,ai,br,bi,cr,ci,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
c=cr+ci*1i;
z=zr+zi*1i;

% Compute relevant Gamma functions using cgama.m [71]
[gr1,gi1]=cgama(ar,ai,1);e1=gr1+gi1*1i;
[gr2,gi2]=cgama(br,bi,1);e2=gr2+gi2*1i;
[gr3,gi3]=cgama(cr,ci,1);e3=gr3+gi3*1i;
[gr4,gi4]=cgama(br-ar,bi-ai,1);e4=gr4+gi4*1i;
[gr5,gi5]=cgama(ar-br,ai-bi,1);e5=gr5+gi5*1i;
[gr6,gi6]=cgama(cr-ar,ci-ai,1);e6=gr6+gi6*1i;
[gr7,gi7]=cgama(cr-br,ci-bi,1);e7=gr7+gi7*1i;
[gr8,gi8]=cgama(cr-ar-br,ci-ai-bi,1);e8=gr8+gi8*1i;
[gr9,gi9]=cgama(ar+br-cr,ai+bi-ci,1);e9=gr9+gi9*1i;

% Compute function using method (a) of Section 4.2 along with (4.16)-(4.20)
if abs(z)<rho
    h=taylora2f1(ar,ai,br,bi,cr,ci,zr,zi,1e-15);
elseif abs(1-z)<rho % Require that |arg(1-z)|<pi
    h=e3*e8/e6/e7...
        *taylora2f1(ar,ai,br,bi,...
        ar+br-cr+1,ai+bi-ci,1-zr,zi,1e-15)...
        +(1-z)^(c-a-b)*e3*e9/e1/e2...
        *taylora2f1(cr-ar,ci-ai,cr-br,ci-bi,...
        cr-ar-br+1,ci-ai-bi,1-zr,zi,1e-15);
elseif abs(z/(z-1))<rho
    h=(1-z)^(-a)...
        *taylora2f1(ar,ai,cr-br,ci-bi,...
        cr,ci,real(z/(z-1)),imag(z/(z-1)),1e-15);
elseif abs(1/z)<rho % Require that |arg(z)|<pi and |arg(1-z)|<pi
    h=(-z)^(-a)*e3*e4/e2/e6...
        *taylora2f1(ar,ai,ar-cr+1,ai-ci,...
        ar-br+1,ai-bi,real(1/z),imag(1/z),1e-15)...
        +(-z)^(-b)*e3*e5/e1/e7...
        *taylora2f1(br-cr+1,bi-ci,br,bi,...
        br-ar+1,bi-ai,real(1/z),imag(1/z),1e-15);
elseif abs(1/(1-z))<rho % Require that |arg(1-z)|<pi
    h=(1-z)^(-a)*e3*e4/e2/e6...
        *taylora2f1(ar,ai,cr-br,ci-bi,...
        ar-br+1,ai-bi,real(1/(1-z)),imag(1/(1-z)),1e-15)...
        +(1-z)^(-b)*e3*e5/e1/e7...
        *taylora2f1(br,bi,cr-ar,ci-ai,...
        br-ar+1,bi-ai,real(1/(1-z)),imag(1/(1-z)),1e-15);
elseif abs(1-1/z)<rho % Require that |arg(z)|<pi and |arg(1-z)|<pi
    h=z^(-a)*e3*e8/e6/e7...
        *taylora2f1(ar,ai,ar-cr+1,ai-ci,...
        ar+br-cr+1,ai+bi-ci,real(1-1/z),imag(1-1/z),1e-15)...
        +z^(a-c)*(1-z)^(c-a-b)*e3*e9/e1/e2...
        *taylora2f1(cr-ar,ci-ai,1-ar,-ai,...
        cr-ar-br+1,ci-ai-bi,real(1-1/z),imag(1-1/z),1e-15);
else
    ['Another method is needed']
    return
end