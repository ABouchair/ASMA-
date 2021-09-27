%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   LDestimate.m                                        %
%                                                       %
%        D. Veitch   P.Abry                             %
%                                                       %
%   1/6/98                                              %
%  DV 4/99                                              %
%  DV 7/99                                              %
%  DV 15/5/2000                                         %
%  DV 8/10/2002 Add extra outputs yj and varj           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function estimates the scaling parameter  'alpha'  of a scaling process.
%    Such processes include  Long Range Dependent (LRD), Self-Similar, 1/f noise, fractal,
%    and multifractal processes.   In the case of LRD a second parameter, cf,  is also estimated,
%    as LRD requires (at least) two parameters for its description.
%      Note:  This function only looks at second order statistics (essentially covariances), this will not give
%             the full picture for multifractal processes.
%
%    More precisely, this function estimates the two parameters of LRD: (alpha,cf), using the
%    wavelet based joint estimator of Abry and Veitch.  It is a spectral domain definition, namely where
%    the spectrum takes the form
%                                          f(nu) ~ cf*nu^alpha near the origin, for a certain range of
%    scales (frequencies).
%
%    In the case where the series to be analysed is intrinsically discrete, a special initialisation procedure
%    must be followed to avoid errors at low scales.  This is activated via one of the arguments.
%
% *** Usage:  [alphaest,cfCest,cfest,Cest,Q,j1opt,yj,varj] = LDestimate(data,regu,j1,j2,discrete_init,calcj1,printout)
%
%--- Routines called Directly:
% wtspec.m
%      " function [muj,nj]=wtspec(data,regu,nbvoies) "
% regrescomp.m
%      " function [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe]
%                                                =  regrescomp(regu, nj, muj,j1,j2,printout) "
% initDWT_discrete.m
%      " function [appro,kfirst,klast] = initDWT_discrete(data,regu,lengthIwant,printout) "
% gauss_CDFinv.m
%      " function [quantile] = gauss_CDFinv(prob, mean, std) "
% newchoosej1.m
%     "  function [j1opt,Qmat] = newchoosej1(regu, nj, muj, printout, j2vec)"
%
%
%  Input:   data:  the input data as a row vector:  this is in fact the sequence of wavelet "approximation"
%                  coefficients, or an approximation thereof (for example often sampled data is used here).
%           regu:  (regularity) =  number of vanishing moments (of the Daubechies wavelet).
%           j1:    the lower limit of the scales chosen,  1<= j1 <= scalemax-1
%           j2:    the upper limit of the octaves chosen, 2<= j2 <= scalemax
%           discrete_init:   1:  perform the special MRA initialisation for intrinsically discrete series
%                         else: assume input data already initialised (ie is already the approximation sequence)
%           calcj1:  1: decision to run the newchoosej1 function, which will output a plot of Q(j_1) 
%                    vs j_1, as well as returning the optimal j1 values (one per j2 input to it).
%                    The values of j2 chosen are not input parameters- need to edit this function.
%                    NOTE: choosing this does Not affect the j1 passed to
%                    LDestimate for estimation. 
%                 else: don't run newchoosej1 
%           printout:   1:
%                        -- a log-log plot graph is plotted plus the regression line and +- 1.96*std(muj)
%                           around each point, being the 95% confidence interval under Gaussian assumptions
%                        -- the values of H etc are printed.
%                        -- a loop is entered allowing interactive choosing of  (j1,j2) (return to exit)
%                    else:   nothing is printed or plotted.
%
%  Output:  alphaest:   estimate of the spectral estimate  alpha
%           cfCest:     estimate of the intermediate quantity  cfC
%           cfest:      estimate of the second LRD parameter  cf = cfC / C
%           Cest:       estimate of the wavelet dependent integral C(alpha) using the estimated alphaest
%           Q:          the goodness of fit measure, the Chi2 based probability that the null hyp is true given
%                       the observed data over the scale range  (j1,j2)
%           j1opt:      the optimal j1^*(j2) vector outputted from newchoosej1  
%           yj:         the LD values
%           varj:       their variances, allowing CI's to be replotted.
%
%   eg of interactive use
%    >> load fgn8.dat                   % discrete data 4096 long, which is  2^12   (12 octaves)
%    >> LDestimate(fgn8,3,1,12,1,1,1);    % try the full range initially, (j1,j2)=(1,12),  print and plot output
%                                       % and initialise as the data is intrinsically discrete
%                                       %  vary (j1,j2) interactively.
%                                       % also run newchoosej1 to get its opinion
%
%   eg of batch use, with (j1,j2)=(4,11),  regu = 2  established by prior experiments.
%   data is sampled real data, so do not performed the discrete initialisation
%    for i= 1:no_realisations
%     [alphaest[i],cfest[i] = LDestimate(realdata(i,:)),2,4,11,0,0,0);     % don't printout,
%                                                     % only store (alpha,cf)
%                                                     % don't initialize,
%                                        % don't calculate optimal values as don't store them
%    end
%--------------------------------------------------------------------------------------------------------
function [alphaest,cfCest,cfest,Cest,Q,j1opt,yj,varj] = LDestimate(data,regu,j1,j2,discrete_init,calcj1,printout)

%%%  Initialize
format compact                          % eliminate excess blank lines in output
n = length(data);    % record original length

%  Massage the input parameters
nbvoies = fix( log2(n) );        % determine the largest possible number of octaves in the data (never realised)

%  Initialize the MRA based recursive method of calculating the DWT coefficients
if  discrete_init==1    % use the special initialisation for intrinsically discrete series
   filterlength = 0;    %  choose the automatic length selection algorithm, or set here if desired
   %  here would be nice to call the output 'appro',  but wastes memory
   [data,kfirst,klast] = initDWT_discrete(data,regu,filterlength,0);          % no output
   if printout
      fprintf('** Using initialization for discrete series, filterlength = %d\n',n+kfirst-klast)
   end
   n = klast-kfirst+1;
else
   if printout
      fprintf('** Taking the given data as the initial approximation sequence\n')
   end
end

%  perform the decomposition using Daubechies wavelets
[muj,nj]=wtspec(data,regu,nbvoies);
if printout
  fprintf('No of points n_j at octave j:   ')
  fprintf('%d ',nj)
  fprintf('\n')
  fprintf('Number predicted by nj=n*2^j:   ')
  %fprintf('%d ',floor(nj(1).*2*2.^(-(1:length(nj)))) )
  fprintf('%d ',floor(n*2.^(-(1:length(nj)))) )
  fprintf('\n\n')
end

 
% Make an automated choice of j1  (value not used in the LD, but separate graph and text displayed)
if (calcj1)
   maxj2 = length(nj);
   %[j1opt] = newchoosej1(regu, nj, muj, printout,4:maxj2);          % full range of possible j2
   %[j1opt] = newchoosej1(regu, nj, muj, printout,[5,8 maxj2]);      % specific choice, optimise this
   [j1opt] = newchoosej1(regu, nj, muj, printout, maxj2);            % LRD choice
   %%% modify plot output of regresscomp
   if discrete_init & printout
      titlehandle = get(gca,'title');
      titlestring = get(titlehandle,'string');   % get current title string
      set(titlehandle,'string',[ titlestring, ',   D-init' ])  % append discrete init indication
   end
else
   j1opt = 1;   % return dummy value to avoid output error
end


% **************** Interactive loop 
while  1
 %  Verify input parameters
 j2 = min(j2,length(muj));    % make sure j2 cannot exceed the maximum number of octaves available
 j2 = max(j2,2);              % make sure that j2>1
 j1 = min(j1,j2-1);           % make sure that j1<j2

 
 %  Perform the joint parameter estimations and calculate the goodness of fit measure, and plot
 [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe,yj,varj] = regrescomp(regu, nj,muj, j1,j2, printout);
 
 %  Calculate the safe octaves for the cfC calculation, and printout a warning if printout=FALSE and some not safe 
 safe_octaves = j1:(j2-unsafe);
 
 if  length(safe_octaves) ~= (j2-j1+1) & ~printout
    fprintf('Warning, %d   octaves  were unsafe, not used to calculate cfC. \n',unsafe)
 end
 
 % Printout a summary of the octaves: in-data/used/available/safe,   and the results on H etc
 if (printout)
  %%% modify plot output of regresscomp
   if discrete_init
      titlehandle = get(gca,'title'); 
      titlestring = get(titlehandle,'string');   % get current title string 
      set(titlehandle,'string',[ titlestring, ',   D-init' ])  % append discrete init indication
   end
   %divergence_test
 
  %%% Text output: 
   fprintf('Octaves:    in_data     available     selected     Goodness of fit (Prob of data assuming lin-regression)\n')
   fprintf('            1--%3.1f       1--%d          %d--%d                      %6.5f \n\n', log2(n), length(nj),j1,j2,Q )
 
  %%%% Calculate Confidence Intervals,   95% two sided. Regresscomp outputs variance, here map to CI
   %  set confidence level.  
   sig_level = 5 ;


   % cf  95% two sided lognormal assumption.  L(M,V)=exp(N(m,v))   m = ln(M^4/(V+M*M))/2;  v = ln(V/M^2 +1)
   z1 = sqrt(2) * erfinv(2*   sig_level/2/100  -1) ;
   z2 = sqrt(2) * erfinv(2*(1-sig_level/2/100) -1) ;
   m = log( (cfest^4)/(Vcf+cfest*cfest) )/2;
   v = log( Vcf/cfest/cfest + 1 );
   z1 = ( z1*sqrt(v) ) + m ;
   z2 = ( z2*sqrt(v) ) + m ;
    cfL = exp(z1);
    cfR = exp(z2);
 
   % Scaling parameters:  95% two sided gaussian assumption 
   seuil=1.9599;
   HLRD = (alphaest + 1)/2;
   H    = (alphaest - 1)/2;
   h    = (alphaest - 1)/2;  
   D    = (5 - alphaest)/2;
   aL = alphaest - seuil*sqrt(Valpha);
   aR = alphaest + seuil*sqrt(Valpha); 
   HLRDL = HLRD - seuil*sqrt(Valpha/4);
   HLRDR = HLRD + seuil*sqrt(Valpha/4);
   HL    = H    - seuil*sqrt(Valpha/4);
   HR    = H    + seuil*sqrt(Valpha/4);
   DL    = D    - seuil*sqrt(Valpha/4);
   DR    = D    + seuil*sqrt(Valpha/4);
  
 
   %  Print the output
   % fprintf('Goodness of fit statistic Q (probability of data assuming regression valid): %6.5f \n \n',Q)
   fprintf('Scaling parameters are:    alpha (LRD)     H (LRD rewrite)  H=h (ss,Holder)  D (frac dim, if alpha in (1,3)) \n')
   fprintf('             Estimates:       %4.3f            (%4.3f)           %4.3f            %4.3f\n',alphaest,HLRD,H,D)
   fprintf('                  CI''s:  [%4.3f, %4.3f]   [%4.3f, %4.3f]   [%4.3f, %4.3f]   [%4.3f, %4.3f]\n\n', aL,aR,HLRDL,HLRDR,HL,HR,DL,DR )

   fprintf(' Second parameters are:         cf              N/A            sigma^2              N/A         \n')
   fprintf('             Estimates:      %7.4f                      Work in progress            \n',cfest)
   fprintf('                  CI''s:  [%7.5f, %7.5f]                Work in progress          \n', cfL,cfR )
   fprintf('\n\n')       

   %  Prompt for new (j1,j2) values
   j1=input('New initial octave j1? (hit return to exit loop)   ');
   if  isempty(j1) 
      hold off
      return
   end
   j2=input('New final octave j2?   ');
   hold off
   fprintf('**************************************************************************************\n');
 else   %  if you can't see the answer there's no point prompting for more values
   return
 end
 hold off
end









