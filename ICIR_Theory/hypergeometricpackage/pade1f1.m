function [h]=pade1f1(ar,ai,br,bi,zr,zi,p,q)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) for parameters%
% a,b and variable z using Pade approximants as detailed in Appendix G.7. %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         p=Number of terms desired on numerator of approximant           %
%         q=Number of terms desired on denominator of approximant         %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Set up vector expressing first p+q+1 terms in Taylor series (3.1)
A=zeros(p+q+1,1);
A(1)=1;
for j1=2:p+q+1
    A(j1)=(a+j1-2)/(b+j1-2)/(j1-1)*A(j1-1);
end

% Set up solution matrix B as in (G.23)
c=zeros(p+q+1,1);
B=zeros(p+q+1,p+q+1);
y=zeros(p+q+1,1);

for j2=1:p+q+1
    c(j2)=A(j2);
end

for j3=1:p+1
    B(j3,j3)=1;
end

for j5=p+2:p+q+1
    for j4=j5-p:p+q+1
        B(j4,j5)=-A(p+1-j5+j4);
    end
end

% Solve matrix system (G.23)
y=B\c;

% Use solution to compile Pade approximant
pvec=zeros(p+1,1);
qvec=zeros(q,1);
numvec=zeros(p+1,1);
domvec=zeros(q,1);

for j6=1:p+1
    pvec(j6)=y(j6);
    numvec(j6)=z^(j6-1);
end

for j7=1:q
    qvec(j7)=y(j7+p+1);
    domvec(j7)=z^(j7);
end

% Return Pade approximant
h=(sum(pvec.*numvec))/(1+sum(qvec.*domvec));