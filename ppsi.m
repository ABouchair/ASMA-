%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   ppsi.m                                               %
%                                                       %
%        D. Veitch   P.Abry                             %
%                                                       %
%   16/06/97                                            %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PPsi = ppsi(x,epsilon)

%%euler= 0.5772156649
%%z=1/2:1:100
%% asymptotic formula for psi(z), z -> inf
%% quite accurate even for z=1 !  not too bad even for z=1/2
%%psi = log(z) - 1/2./z - 1/12./z./z + 1/120./z.^4 - 1/252./z.^6
%
%   x is a matrix (or vector, or single value) of values at which the 
%       Psi or DiGamma function  psi(x) is to be calculated.
%   epsilon is some kind of absolute  error.  
%   
%   If x<N then x is first increased up to x+N using the  recursion relation 
%      Psi(z) = Psi(1+z) - 1/z, 
%   and then the much LARger value of Psi(x+N) is  calculated  with an 
%   asymptotic formula.  N is determined by epsilon and a guessed order
%   of magnitude of the error in the asymptotic formula. 
%
% file:      	psi.m
%

N = ceil( 252*epsilon^(-1/6) );

qs =size(x);
for i=1:qs(1)
  for j=1:qs(2)
    z= x(i,j);
    if (z>N)
       PPsi(i,j) =  log(z)     - 1/2/z     - 1/12/z^2     + 1/120/z^4     - 1/252/z^6;
    else
       PPsi(i,j) =  log((z+N)) - 1/2/(z+N) - 1/12/(z+N)^2 + 1/120/(z+N)^4 - 1/252/(z+N)^6            - sum( (z+(0:N-1)).^(-1) );
    end
  %fprintf(1,'%12.5f  %8d %24.10e   %24.10e \n',x(i,j), N, log(z)     - 1/2/z     - 1/12/z^2     + 1/120/z^4     - 1/252/z^6, Psi(i,j));
  end
end
