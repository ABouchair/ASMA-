function pquantile = chi2_CDFinv(p,N)
%
%   p est la probabilite pour la laquelle on cherche la quantile correspondante.
%   N no de degre de liberte du loi de Chi2.
%
%  pour Chi2_N, le CDF est   p = F(x) =  gammainc(x/2,N/2)
%               On cherche   x = F^-1(p), le p'ieme quantile.
%
%   on passe  F(x) - p  `a  fzero

% cherche pour le zero dans l'interval  x\in [0,20*N]
% pquantile = fzero('chi2_CDFmod',[0,20*N],[],[],N,p);  % l'ancienne version de Matlab
pquantile = fzero('chi2_CDFmod', [ 0 20*N], optimset([]), N,p);  %  

%  Merde cette Subfunction pour CDF marche mais fzero ne le reconnais pas.
%function prob = CDF(x,N,p)
%prob = gammainc(x/2,N/2) ;    % vrai CDF
%prob= prob - p;               % modifie pour fzero

