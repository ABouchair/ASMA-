%  This function inverts the Gaussian PDF to yield the
%  p-quantile of a normal rv with mean m and variance sig^2.
%
%      quantile = gauss_CDFinv(p,m,sig)
%
function  quantile = gauss_CDFinv(p,m,sig)

quantile = sig*(  sqrt(2) * erfinv(2*p-1) )  +  m  ;
