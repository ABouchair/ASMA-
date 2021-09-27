%  This function gives values from  the Gaussian PDF
%  with mean m and variance sig^2.
%
%  A vector of means for constant variance will work
%  or a vector of means and variances.
%  Vector q does not work however
%
%      prob = gauss_CDF(q,m,sig)
%
function  prob = gauss_CDF(q,m,sig)

%fprintf('%f %f\n',m,sig)
prob = ( erf(  (q-m)./sig   /sqrt(2)) + 1 )/2  ; 
