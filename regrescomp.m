%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%  regrescomp.m                                         %
%                                                       %
%       P. Abry and D. Veitch                           %
%                                                       %
%   17/07/97                                            %
% DV, Melb  15 May 2000                                 %
% DV, Melb, 10 Sep 2004:  robustified j1,j2 choice      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%--- Routines used directly:
% rzeta.m
% psi.m
% rsynthond1.m
% ffteur.m
%
% regrescomp.m estimates alpha, cfC, C and cf using the weighted linear regression with the bias "modification".
%          It fully implements the Joint estimator paper, and calculates
%          confidence intervals and variances from expressions, not from data
%
%  Input:   regu:  number of vanishing moments (of the Daubechies wavelet).
%           nj:   the vector[1..scalemax] of the actual number of coefficients averaged at octave j. (these
%                  numbers are not powers of two, even if the length of the data is, because of the border
%                  effects at each scale). 
%           muj:   the vector[1..scalemax] of the average of the squared wavelet coefficients
%     scalemax  is NOT necessarily the largest octave such that there are still some coefficients left. In addition
%             extra octaves are removed on the right before passing them to this routine, if the number of coefficients is
%             equivalent to less than twice the "length" of the wavelet used.
%               
%            j1:   the lower limit of the scales chosen,  1<= j1 <= scalemax-1  
%            j2:   the upper limit of the octaves chosen, 2<= j2 <= scalemax
%            printout:  if = 1  then printout a log-log plot plus the regression line and +- 2*std(muj) around each point
%           
%  Output:  alphaest:   estimate of the spectral parameter  alpha
%           cfCest:     estimate of the intermediate quantity  cfC
%           cfest:      estimate of the second LRD parameter  cf = cfC / C
%           Cest:       estimate of the wavelet dependent integral C(alpha) using the estimated alphaest
%           Q:          the goodness of fit measure, the Chi2 based probability that the null hyp is true given the
%                       observed data over the scale range  (j1,j2)
%           vjj:       the coefficients for the intercept of the weighted linear fit, for use in the testing of
%                       the safety of the scales used in the cfC calculation 
%
%
%        Covariance matrices for the weighted estimators, under the hypotheses of Gaussianity and complete decorrelation.
%           Valpha:        Exact variance of alphaest  
%           VcfC:          Estimator of exact variance of cfC   (biased, ~ square of cfC)
%           CoValphacfC:   Estimator of exact covariance of (alphaest,cfCest)  (unbiased! linear fn of cfC)
%           Vcf:           Estimator of approximate variance of cf                   
%           CoValphacf:    Estimator of approximate covariance of (alphaest,cfest) 
%
%        unsafe:       the number of scales which were found to not pass the test for the existence of the variance of cfC.
%           yj:            The dependent variables of the regression
%           varj:          The variances of yj
%           aest:          The intercept of the linear regression, a-estimate 
% 
%   *** Note,  j   ranges over 1..scalemax    nj, muj, yj, varj  etc
%              jj  ranges over j1..j2         njj,mjj,yjj,varjj  etc
%
%  function [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe,yj,varj,aest]=regrescomp(regu, nj, muj,j1,j2,printout);
%
%-----------

function [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe,yj,varj,aest]=regrescomp(regu, nj, muj,j1,j2,printout);


%%  Internal parameters
fsize = 14;      % set font size (20 for papers)

%--- Selection confidence interval type 
ICtype = 1;        % Gaussian only
ICtype = 2;        % Data based only
ICtype = 3;        % Both

%--- Constants 
precis_zeta = 1e-4 ;
precis_psi  = 1e-4 ;
loge  = log2(exp(1)) ;
 
%--- Variables
scalemax = length(muj);
j2 = min(j2,scalemax);   % make sure j2 is not chosen too large
if (j1>=j2) 
   fprintf(' ** Regrescomp:  scale range given: (%d,%d) is only 1 wide! impossible to estimate\n',j1,j2);
   alphaest=0; cfCest=0; cfest=0; Cest=0; Q=0;
   Valpha=0; VcfC=0; CoValphacfC=0; Vcf=0; CoValphacf=0; unsafe=0; yj=0; varj=0; aest=0;
   return;
end
jj = j1:1:j2;
J = length(jj); 
njj  =  nj(jj);
mujj = muj(jj);

%--- Compute the bias correction for all octaves
 gj = ppsi(nj/2,precis_psi) * loge - log2(nj/2) ;
 
% -- Calculate the random variables of the linear regression for all
% octaves, then for  j1..j2
%  muj
%  gj
yj = log2(muj) - gj ;
yjj = yj(jj) ;
 
%--- Compute la variance des yj = log(muj) - g(j)  for all octaves, then for  j1..j2
varj = loge^2 * rzeta(fix(nj/2),precis_zeta) ;  
varjj = varj(jj) ;


%--- Fourier transform of the wavelet  (used in the calculation of Cest)
nbvoies3=6;
wsf=2;
nbpoints2=16384;
[approx,tapprox,napprox]=rsynthond1(regu,wsf,nbvoies3);
[f1,f2]=ffteur(approx,nbpoints2,2^nbvoies3);
fe=f2(2)-f2(1);
fr=f2(1:nbpoints2/2);fr=-fliplr(fr);
ft=f1(1:nbpoints2/2);ft= fliplr(ft);

%--- Compute the coefficients for octaves  j1..j2
S0 = sum(1./varjj) ;
S1 = sum(jj./varjj) ;
S2 = sum(jj.^2./varjj) ;
wjj = (S0 * jj - S1) ./ varjj / (S0*S2-S1*S1) ; 
vjj = (S2 - S1 * jj) ./ varjj / (S0*S2-S1*S1) ;

%%  Test la calculation des vjj et des wjj
%sum(jj.*wjj)
%sum(vjj)
%sum(wjj)
%sum(jj.*vjj)

%   Select those octaves which are safe for the cfC calculation (those where the variance exists)
%   Careful, these indices should only be applied to vectors which are already of the form jj = j1..j2
jj_safe  = find( (njj+4*vjj)  > 0.02);      
if length(jj)~=length(jj_safe)   
    fprintf('\n*** Regrescomp: range choice (j1,j2)=(%d,%d) not safe for cfC and cf, regression will have to be redone!).\n',j1,j2);  
    %fprintf('Estimation for alpha is reliable however.   Test on octaves gives: nj+4*vjj = \n\n');
    njj+4*vjj;
    unsafe = length(jj) - length(jj_safe);
else
    unsafe=0;
end    
   

%  diagnostics
%jj
%jj_safe
%vjj
%find( (njj+4*vjj) > 0.02)

%  estimate  a, b=alpha, and C 
alphaest  = sum(wjj .* yjj ) ;
aest      = sum(vjj .* yjj ) ;
cfCest    = prod(mujj.^vjj ) ;
Calpha = alphaest;
if (alphaest<0)
  %fprintf('**  alpha negative with (j1,j2)=(%d,%d), set to 0 for Cest calculation \n',j1,j2)
  Calpha=0; 
end
if (alphaest>1)
  %fprintf('**  alpha > 1 with (j1,j2)=(%d,%d), set to 1 for Cest calculation \n',j1,j2) 
  Calpha=1;
end
Cest      = 2*sum(fr.^(-Calpha)  .* (abs(ft)).^2)*fe ;        % use the wavelet shape calculation here

%--- corrective multiplicative bias factor for cfC
%pmu = exp( sum( gammaln(njj(jj_safe)/2)  +  vjj(jj_safe).*log(njj(jj_safe)/2)  -  gammaln(vjj(jj_safe) + njj(jj_safe)/2)  )) ;
p = exp( sum( gammaln(njj(jj_safe)/2) + vjj(jj_safe).*log(njj(jj_safe)./(2.^(1-gj(jj_safe)))) - gammaln(vjj(jj_safe) + njj(jj_safe)/2) ));
pmu = p * 2^(-sum(vjj.*gj(jj)));

%  estimate  cfC and  cf.
cfCest  = cfCest*pmu;
cfest   = cfCest  / Cest ;


%  Calculation of the theoretical Cov function of (alphahat,cfChat), and the estimate for (alphahat,cfhat)
%  using the estimated values of the (alpha,cfC,cf) 

% ####   alpha
Valpha = sum(varjj.*wjj.*wjj);

% #####  cfC
VcfC = cfCest^2*(exp( sum( gammaln(2*vjj(jj_safe)+njj(jj_safe)/2)+gammaln(njj(jj_safe)/2)-2*gammaln(vjj(jj_safe)+njj(jj_safe)/2) ))-1);
CoValphacfC  =  cfCest * loge * sum(wjj(jj_safe) .* ( ppsi(njj(jj_safe)/2+vjj(jj_safe),precis_psi)-ppsi(njj(jj_safe)/2,precis_psi)) ) ;


% #####  cf
%  cf1est = cfCest/C1(alphaest) = g(x,y) = x*h(y)  ( here x== cfC, y==alpha )
%  Calculate partial derivatives of g(x,y), evaluated at (cfC,alpha) [in fact at (cfCest,alphaest) ]
%  g'=[h,xh'],   g''=[0,   h' ]
%                    [h'  xh'']
pal = 2^alphaest;
h   = ( 1 - alphaest )/( 2 - pal ) ;
hd  = ( log(2)*pal*h -1 )/( 2 - pal ) ;
hdd = log(2)*pal*( h*log(2) + 2*hd )/( 2 - pal ) ;
gx = h ;
gy = cfCest*hd  ;
gxy = hd ;
gyy = cfCest*hdd  ;
%gx*gx*VcfC 
%2*gx*gy*CoValphacfC 
%gy*gy*Valpha 

% ***** Examine the size of the correction factor for  E(cfhat)
%cf_biais = ( 2*gxy*CoValphacfC + gyy*Valpha )/2
%cfest   = cfest + cf_biais;
%cf_biais_relative = cf_biais/cfest

Vcf = gx*gx*VcfC + 2*gx*gy*CoValphacfC + gy*gy*Valpha ;
CoValphacf = h*CoValphacfC + cfCest*hd*Valpha;


%%%%--- Goodness of fit
%   If at least 3 points, Apply a Chi-2 test , no level chosen, should be viewed as a function of j1

if J>2
  J2 = (J-2)/2 ;
  X  = sum( ((yjj - alphaest  * jj - aest ).^2)./ varjj ) ;
  Q  = 1  - gammainc(X/2,J2) ;

 %  alternate calculation taking the differences first then using the test statistic of A4 
  if  0 
    XX = diff(yjj);   % series of samples with the same mean (under H0)
    for i = 1:J-1      % assume the yj are independent so can add variances
      varXX(i) = varjj(i) + varjj(i+1);
    end
    S0XX = sum(1./varXX) ;
    S1XX = sum(XX./varXX) ;
    V = sum( ((XX - S1XX/S0XX).^2)./varXX );
    Qtest = 1 -  gammainc(V/2, J2)
  end
else
  fprintf('\n***** Regrescomp: Cannot calculate Q, need at least 3 points.\n');
  Q=0;
end
  

%%    Plot des resultats if argument "printout" is true
if (printout)
   %fprintf('n/2=nj(1)= %d, whereas the product of (nj^vj)/2 = %d \n', nj(1), floor(prod( (njj/2).^vjj )))
   %fprintf('The corrective bias factor p = %6.3f \n\n',p)
   
   % Calculate confidence intervals
   seuil=1.9599;   % .975 quantile of standard Gaussian
   CI(1,:) = yj-seuil*sqrt(varj);   % Gaussian CI's
   CI(2,:) = yj+seuil*sqrt(varj);
   % Determine extremes for plotting
   mi = min(CI(1,:)) - 0.2;  
   ma = max(CI(2,:)) + 0.2;  

   figure(11);  clf
   plot(yj)
   %axis([j1-0.5 j2+0.5 mi ma])     % focus in on the octaves chosen
   axis([0.8 scalemax+0.2 mi ma])     % show the confidence intervals to either side 
   hold
   set(gca,'FontSize',fsize-2)
   %plot(yj,'*')

   %  Plot the vertical bars at each octave showing seuil* the standard deviation of the y(j)
   for k=1:1:scalemax
      plot([k k ], CI(:,k) ) ;
   title(['Logscale Diagram,  N=',num2str(regu),'    [ (j_1,j_2)= (',num2str(j1),',', num2str(j2),        '),   \alpha-est = ',num2str(alphaest,3), ',    Q= ',num2str(Q),' ]'])
   end
  
   plot(jj,alphaest * jj + aest,'r')
   xlabel('Octave j','FontSize',fsize)
   %ylabel('y(j) = log2( muj ) - g(j) )')
   ylabel('y_j ','Rotation',0,'FontSize',fsize)   % make it upright (doesn't always work!
   grid on; hold off

end





