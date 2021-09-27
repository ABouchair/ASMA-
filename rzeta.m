%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% file:      	rzeta.m
% created: 	Fri May 23 1997 
% authors:  	Matthew Roughan  Darryl Veitch
% email:   	matt@serc.rmit.edu.au
%
%
%   q is a matrix (or vector, or single value) of values at which the
%       generalised Riemann Zeta function  Z(2,q) is to be calculated.
%   epsilon is an upper bound on the RElative error. 10^-6 is around the
%       machine precision. 10^-5 is excellent.
%   For large q,   Z(2,q) ~ 1/q
%   For small q, the first N(q,epsilon) terms are calculated explicitly, then the
%   tail is estimated with relative precision of epsilon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zeta = rzeta(q,epsilon)

N = max(round(epsilon^(-1) -q -1) +1, 0);

qs =size(q);
for i=1:qs(1)
  for j=1:qs(2)
    zeta(i,j) = sum( (q(i,j)+(0:N(i,j))).^(-2) ) + (q(i,j)+N(i,j)+1)^(-1);
    %fprintf(1,'%12.0f  %8d %24.10e   %24.10e \n',q(i,j), N(i,j), 1/q(i,j), zeta(i,j));
  end
end
