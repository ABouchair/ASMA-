% version modifiee pour aider la function chi2_CDFinv
function  prob = chi2_CDFmod(x,N,p)
prob = gammainc(x/2,N/2) ;    % vrai CDF pour Chi2_N(x)
prob= prob - p;               % modifie pour fzero
