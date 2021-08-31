function [h]=bessel0f1(ar,ai,zr,zi)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric limit function 0F1(;a;z) using    %
% the Bessel functio representation stated in Appendix I.                 %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
% Output: h=Computed value of 0F1(;a;z)                                   %
%-------------------------------------------------------------------------%

% Compute a,z in terms of ar,ai,zr,zi
a=ar+ai*1i;
z=zr+zi*1i;

% Compute Bessel function representation (note that due to the lack of
% capabilities of 'besselj.m', this can only work if 2i*sqrt(z) is real)
h=gamma(a)*z^((1-a)/2)*1i^(1-a)*besselj(a-1,2i*sqrt(z));