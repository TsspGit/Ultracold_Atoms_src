function [h]=singlefraction1f1(ar,ai,br,bi,zr,zi,tol)

%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the confluent hypergeometric function 1F1(a;b;z) using the     %
% single fraction method of Section 3.3.                                  %

%-------------------------------------------------------------------------%
% Input:  ar=Re(a)                                                        %
%         ai=Im(a)                                                        %
%         br=Re(b)                                                        %
%         bi=Im(b)                                                        %
%         zr=Re(z)                                                        %
%         zi=Im(z)                                                        %
%         tol=Specified tolerance as detailed in Section 3.3              %
% Output: h=Computed value of 1F1(a;b;z)                                  %
%-------------------------------------------------------------------------%

% Compute a,b,z in terms of ar,ai,br,bi,zr,zi
a=ar+ai*1i;
b=br+bi*1i;
z=zr+zi*1i;

% Start with 2-vector of values a1(1),a1(2),b1(1),b2(2),c1(1),c1(2),d1(1),
% d1(2), denoted as alpha(0),alpha(1),beta(0),beta(1),gamma(0),gamma(1),
% delta(0),delta(1) in Section 3.3
a1=[0,b];b1=[1,a*z];c1=[1,b];d1=[1,(b+a*z)/b];

for j=3:500
    % Compute next values of a1(j),b1(j),c1(j) in terms of last ones
    a1(j)=(a1(j-1)+b1(j-1))*(j-1)*(b+j-2);
    b1(j)=b1(j-1)*(a+j-2)*z;
    c1(j)=c1(j-1)*(j-1)*(b+j-2);
    % Stop if any of these are infinitely large, so MATLAB won't be able to
    % carry out addition or division with them
    if (a1(j)==Inf)||(b1(j)==Inf)||(c1(j)==Inf)
        break
    end
    % Compute next term in d1(j)
    d1(j)=(a1(j)+b1(j))/c1(j);
    % If stopping criterion as detailed in Section 3.3 has been satisfied,
    % terminate computations
    if abs(d1(j)-d1(j-1))/abs(d1(j-1))<tol && abs(d1(j-1)-d1(j-2))/abs(d1(j-2))<tol
        break
    end
    % If 500 terms have been computed without stopping criterion being
    % satisfied, state this
    if (j==500)
        [' ' num2str(j) ' terms computed']
        return
    end
end

% If stopping criterion has been satisfied, return sum of terms computed,
% i.e return d1(end)
h=d1(end);