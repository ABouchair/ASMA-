%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   eda_staty.m                                         %
%                                                       %              
%        D. Veitch   P.Abry                             %
%                                                       %
%   25/1/98                                             %
%   DV, Lyon  14 oct                                    %
%   DV, Melbourne, 2 May 2000                           %
%   DV, Melbourne, 24 Dec 2000                          %
%   DV, Melb, 10 Sep 2004:  robustified j1,j2 choice    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function supports an experimental analysis of stationarity. 
%    The series is cut into subseries and tests performed on each for:
%        - the mean   , with confidence intervals: LRD and IID
%        - the variance
%        - H   and the stationarity test  (alpha or H are possible, see the internal variable "wantH" 
%                                          which selects which one)
%        - cf
%    
%    Other choices can easily be added!! 
%
%             eg:  eda_staty(delay45,1,1,1,1,1,20,  2   ,2,6,2,6, 1,1);   
%                  Examine delay45  cut into twenty equal size 
%                  subseries, considering mean, var, and H with N=2, and plotting the series
%                  Test is for the regime [j1,j2]=[2,6]
%                  eda_staty(fgn8     ,1,1,0,1,1, 5,  2   ,2,7,2,7  , 1, 1);  % with initialisation j1*=2 !
%                  eda_staty(data     ,1,1,0,0,0,  6,  2   ,6,13,6,13, 0, 1);  % demo for CI's of mean
%                  eda_staty(pOct_work',0,1,0,1,0, 4,  2   ,6,14,6,14, 0, 1);  % acceptance of constancy of alpha     
%
%
%--- Routines called Directly :
%
%   % wtspec, newchoosej1,  regrescomp, initDWT_discrete, chi2_CDF, chi2_CDFinv, gauss_CDF, gauss_CDFinv
%
%  Input:   x:    the series to be analysed (is converted to a row vector)
%           looksig :   logical parameter, 1 means the signal will be plotted
%           lookmean:   logical parameter, 1 means the mean will be  analysed, else should be zero.
%           lookvar :   logical parameter, 1 means the variance   "   "
%           lookH   :   logical parameter, 1 means the exponent (H or alpha)   "   "
%           lookcf  :   logical parameter, 1 means the cf S       "   " 
%           cutpoints:  value or vector describing how the series should be partitioned.
%                         - if a value=nsub, then nsub subseries are created of equal size.
%                         - else a vector of specific values is taken (*** not currently supported***) 
%           regu:  (regularity) =  number of vanishing moments (of the Daubechies wavelet).
%           j1slice:   the j1 value for each block
%           j2slice:   the j2 value for each block
%           j1full:   the j1 value for the entire series
%           j2full:   the j2 value for the entire series
%           discrete_init: 1:  perform the special MRA initialisation for intrinsically discrete series
%                         else: assume data already initialised (ie is already the approximation sequence)
%           printout:   1: -text output and the multiplot graph, one row per statistic selected.
%                          -if LD's are calculated, then outputted in a plot, one for each block, Q in title.
%                    else: nothing.
%
%  Output:  m:   the vector of means of each subseries
%           v:   the vector of variances of each subseries
%           al:   the vector of scaling exponents (H or alpha) values of each subseries 
%           cf:    the vector of cf values of each subseries
%           stdal: the corresponding vector of standard deviations of the estimates (theoretical)
%           stdcf: the corresponding vector of standard deviations of the estimates (theoretical)
%           fullal:     the exponent value (H or alpha) over the entire series
%           fullcf:     the cf value over the entire series
%           fullstdal:  the corresponding  standard deviation of the estimate (theoretical)
%           fullstdcf: the corresponding  standard deviation of the estimate (theoretical)
%           crit_level:  the 'complementary' quantile (expressed as a probability) of the constantcy test statistic, ie probability that the data is at least this bad.
%           crit_level2: a vector of values for the two-level test with m1 = 1..nsub/2  .
%
%  A set of subplots are plotted on figure 20.  First the series, if selected, followed by 
%  a subplot for each quantity selected 
%
%  Call: [m,v,al,cf,stdal,stdcf,fullal,fullcf,fullstdcf,fullstdal,crit_level,crit_level2] = eda_staty(x, looksig, lookmean, lookvar, lookH, lookcf, cutpoints , regu, j1slice, j2slice, j1full, j2full, discrete_init, printout)
%

function [m,v,al,cf,stdal,stdcf,fullal,fullcf,fullstdcf,fullstdal,crit_level,crit_level2] = eda_staty(x, looksig, lookmean, lookvar, lookH, lookcf, cutpoints , regu, j1slice, j2slice, j1full, j2full, discrete_init, printout)

%%%%%%%%%% internal variables
wantH = 0;   % print and output according to H or alpha,   alpha = 2H - 1


%%%%%%%%%%% correct input
if  (lookH | lookmean | lookcf)  % if want to look at H, make sure regu is set.
   regu=max(regu,1);
end
dim = size(x);
if dim(1) > 1
  x = x';
end

%%%%%%%%%%% init
x = x';
len = length(x);
m = 0;
v = 0;
al = 0;
stdal = 0;
fullal = 0;
fullstdal = 0;
crit_level = 0;


%%%%%%%%%%%  determine subseries of x
if (length(cutpoints)==1)
  %cutpoints = len*(1:cutpoints)/cutpoints
  %nsub = length(cutpoints)
  nsub  = cutpoints;
  lensub = floor(len/nsub);
  %cutpoints
  %endpts =  (1:nsub)*lensub ;    % output indices of endpoints of each block
  %endpts(endpts==4424)  % check special value for A4
else
  nsub = length(cutpoints) + 1;
  lensub = diff([1 cutpoints len]) % *** after this unequal lengths are not supported! but could be done
end 
   

%%%%%%%%%%  create subseries as rows of matrix 'series'
series = ones(nsub,lensub);
for n = 1:nsub
  series(n,:) = x( (1+(n-1)*lensub):(n*lensub) )';
end

series(:,1:min(12,lensub));

if ( (lookmean | lookvar | lookH | lookcf ) & printout)
fprintf('\n*******************************************************************************************\n\n')
  fprintf('%d subseries of length %d will be examined\n\n',nsub,lensub)
end


%%%%%%%%%%%  For each quantity, calculate the value for each subseries
seuil = 1.9599 ;
num_plot = ( looksig + lookmean + lookvar + lookH + lookcf );
nplotnum = (num_plot)*100 + 10;
if (printout)
  figure(20)
  clf
end
if (looksig & printout)
   nplotnum = nplotnum+1;
   subplot(nplotnum)
   skip = ceil(len/5000);
   plot(1:skip:len,x(1:skip:len),'k-');   % don't plot all the points of huge timeseries
   axis tight
   %title('Timeseries:  fGn of length 2^{18},  \alpha = 0.6,  c_f = 100')
   title('Timeseries')
   %title('Timeseries, pOct aggregated in 10ms bins')
end

line   = ones(1,nsub);     
CIline = ones(1,nsub+1);

if (lookH | lookcf | lookmean)   %  need H values For confidence intervals for mean as well as for H itself
   alpha = [];
   H = [];
   cf = [];
   subj1 = [];
   subQ = [];

   if (printout)
     fprintf('\n\n');
     figure(21)
     clf
     minlower = 100000000;    %  variables to fit vertical axis in LD plots
     minlowerfull = minlower;
     maxupper= -minlower;
     maxupperfull = maxupper;
   end

   for n = 1:nsub
      %  Initialize the MRA based recursive method of calculating the DWT coefficients
      if  discrete_init==1    % use the special initialisation for intrinsically discrete series
         filterlength = 0;    %  choose the automatic length selection algorithm, or set here if desired
         [appro,kfirst,klast] = initDWT_discrete(series(n,:),regu,filterlength,0);      % no output
         if printout
            fprintf('** Using initialization for discrete series, filterlength = %d\n',lensub+kfirst-klast)
         end
         %n = klast-kfirst+1;
         [muj,nj]=wtspec(appro,regu,fix( log2(klast-kfirst+1) ));    % use length of appro, after filtering
         if printout 
           figure(21);
         end 
      else
         appro = series(n,:);
         if printout
            fprintf('** Taking the given data as the initial approximation sequence\n')
         end
         [muj,nj]=wtspec(series(n,:),regu,fix( log2(lensub) ));
      end

   %  j2 = length(nj);
      j2 = min(j2slice,length(nj)) ;
      if (j2<2 )
        fprintf('Subseries too short to evaluate H! aborting.  \n');
        %j1 = j1slice;
        %H = 1/2.*line;                                                                  
        return;
      else
%        j1 =  newchoosej1(regu,nj, muj, 0,j2 );  % no printout, just the LRD case.
        if (j1slice>=j2)
          fprintf('Subseries too short to evaluate H! please change (j1,j2) selection. For now, trying j1-1 \n');
          j1 = max(j2-1,1);
        else
          j1 = j1slice;
        end
        subj1 = [ subj1 j1 ];
        [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe,yj,varj,aest]= regrescomp(regu,nj,muj, j1,j2,0);
        alpha = [ alpha alphaest ];
        H  = [ H (alphaest + 1)/2 ];
        cf = [ cf cfest ];
        subQ = [ subQ Q ];
        stdalpha(n) = sqrt(Valpha) ;
        stdH(n) = stdalpha(n)/2;
        stdcf(n) = sqrt(Vcf);

        if (printout)
          if (wantH)
            fprintf('For subseries %d, (j1,j2)= (%d,%d), H = %4.2f, cf = %7.4f, Q=%5.4f .\n',n,j1,j2,H(n), cf(n),Q)
          else
            fprintf('For subseries %d, (j1,j2)= (%d,%d), alpha= %4.2f, cf = %7.4f, Q=%5.4f .\n',n,j1,j2,alpha(n),cf(n),Q)
          end 
          %%%  plot the next logscale diagram (in a row.)
          subplot(1,nsub,n)
          plot(yj,'*-')
          hold
          grid
          title(['Q = ',num2str(Q)])
          %title(['\alpha = ',num2str(alphaest,2),',  Q = ',num2str(Q)])
          %xlabel('Octave j')
          %ylabel('y(j)')
          upper = yj+seuil*sqrt(varj);
          lower = yj-seuil*sqrt(varj);
          jj = j1:j2;
          maxupperfull = max(maxupperfull,max(upper));
          minlowerfull = min(minlowerfull,min(lower));
          maxupper = max(maxupper,max(upper(jj)));
          minlower = min(minlower,min(lower(jj)));
          for k=jj    % plot vertical confidence intervals
             plot([k k ],[lower(k) upper(k)]) ;
          end
          plot(jj,alphaest * jj + aest,'r')
          %V = axis;
          %axis([0.8 j2+0.2  V(3) V(4) ])
          hold off
        end
      end
   end
   %  set LD's to the same vertical scale:
   if printout 
      for n = 1:nsub
        subplot(1,nsub,n)
        axis([j1-0.2 j2+0.2 minlower maxupper]);
        %axis([0.2 length(nj)+0.6 minlowerfull-0.5 maxupperfull-1.5]);% special for slide: full j range
      end
   end
   CIalpha = alpha;
   CIalpha(CIalpha<0)=0.01; %  if alpha -ve set to 0
   CIalpha(CIalpha>1)=.9;   %  if alpha >1  set to 1   
   cgam = cf./( 2*((2*pi).^CIalpha) .* gamma(CIalpha) .* sin(pi*(1-CIalpha)/2) );    %  need this for IC's for mean
   LRDasympstd = sqrt(cgam)./( lensub.^((1-CIalpha)/2) .* sqrt((1+CIalpha).*CIalpha/2) );  % asymptotic LRD CI's
   
   %%% Calculate and store H and cf over the entire series,  and the corresponding asympstd for the mean CI.
   %  Initialize the MRA based recursive method of calculating the DWT coefficients, then get the details 
   if  discrete_init==1    % use the special initialisation for intrinsically discrete series
      filterlength = 0;    %  choose the automatic length selection algorithm, or set here if desired
      [appro,kfirst,klast] = initDWT_discrete(x,regu,filterlength,0);      % no output
      if printout
         fprintf('** Using initialization on the full series, filterlength = %d\n',len+kfirst-klast)
      end
      %n = klast-kfirst+1;
      [muj,nj]=wtspec(appro,regu,fix( log2(klast-kfirst+1) ));
      clear appro                                       % this can be big! 
   else
      if printout
         fprintf('** Taking the given data as the initial approximation sequence\n')
      end
      [muj,nj]=wtspec(x,regu,fix( log2(len) ));
   end

   %j2 = length(nj);
   %j1 =  newchoosej1(regu,nj, muj, 0,j2 );
   %j1= max(subj1);
   j1=j1full;
   j2=j2full;
   [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe,yj,varj,aest] = regrescomp(regu,nj, muj, j1,j2,0);
   fullalpha = alphaest;
   fullstdalpha = sqrt(Valpha);
   fullH = (alphaest + 1)/2;
   fullstdH = sqrt(Valpha)/2;  
   fullcf = cfest;
   fullstdcf = sqrt(Vcf);
   CIalpha = fullalpha;
   CIalpha = max(CIalpha,0.01);
   CIalpha = min(CIalpha,0.9);
   cgam = cfest./( 2*((2*pi).^CIalpha) .* gamma(CIalpha) .* sin(pi*(1-CIalpha)/2) )   ;
   LRDasympstdfull = sqrt(cgam)./( lensub.^((1-CIalpha)/2) .* sqrt((1+CIalpha).*CIalpha/2) )  ;
   
   if wantH   % set output variables correctly
      al = H;
      stdal = stdH;
      fullal = fullH;
      fullstdal = fullstdH;
   else
      al = alpha;
      stdal = stdalpha;
      fullal = fullalpha;
      fullstdal = fullstdalpha;
   end 

   if (printout)
     figure(20)
   end
end
   


if (lookmean)  
   m = mean(series');
   stdm = std(series');
% output for Balaton,
seuil*stdm(1)/sqrt(lensub)
seuil*LRDasympstd(1)
   if (printout)
     nplotnum = nplotnum+1;
     subplot(nplotnum)
     plot(1:nsub,m,'k-')
     hold
     plot(1:nsub,mean(m).*line,'k-')
     for n = 1:nsub
%seuil*stdm(n)/sqrt(lensub)
       plot([n n],[m(n)-seuil*stdm(n)/sqrt(lensub) m(n)+seuil*stdm(n)/sqrt(lensub) ],'b-') ;   % add SRD confidence intervals
       if ( H(n)> 1/2 & H(n) <1)
         plot([n+.05 n+.05],[m(n)-seuil*LRDasympstd(n) m(n)+seuil*LRDasympstd(n) ],'r-');% add asymptotic LRD CI's
       else
         plot([n+.05 n+.05],[m(n)-seuil*stdm(n)/sqrt(lensub) m(n)+seuil*stdm(n)/sqrt(lensub) ],'r-');  % SRD CI's if Hest<1/2.
       end
     end
     V = axis;
     axis([ 0.5 nsub+0.5 V(3) V(4) ]); 
     %xlabel('Subseries number')
     %xlabel('Block number')
     %ylabel('Mean rates, megabits/s')
     ylabel('Means')
    % title('Means of the subseries. The horizontal lines give the overall value and confidence intervals.')
     title('Means over the blocks. The horizontal line gives the overall mean.')
     grid

     %  Add CI's for entire trace on the left of mean plot.
     meanfull = mean(m) ;
     stdfull  = std(x)  ;
     plot([.7 .7],[meanfull-seuil*stdfull/sqrt(len) meanfull+seuil*stdfull/sqrt(len) ],'b-');  % SRD to compare
     if ( fullH> 1/2 & fullH <1)
       plot([.75 .75],[meanfull-seuil*LRDasympstdfull meanfull+seuil*LRDasympstdfull ],'r-');   % LRD
     %  plot([.8,1:nsub],meanfull+seuil*LRDasympstdfull.*CIline,'r-')  %  top line of CI 
     %  plot([.8,1:nsub],meanfull-seuil*LRDasympstdfull.*CIline,'r-')  %  bottom line of CI
     else
       plot([.7 .7],[meanfull-seuil*stdfull/sqrt(len) meanfull+seuil*stdfull/sqrt(len) ],'b-');% SRD if H<1/2 or>1
     end
     hold off
   end
end
if (lookvar)
   v = (std(series')).^2;
   vfull = std(x)^2;
   %  if X~N(mu,v), then the sample var (biased) is ~ v/n * ChiSq_n-1, with var = 2v^2 *(n-1)/n^2
   %  assume then that for large n the sample var is Gaussian with var 2v^2/n
   stdasympSRD = sqrt(2/lensub)*v ;   % using v to estimate the real v for each subseries
   if (printout)
     nplotnum = nplotnum+1;
     subplot(nplotnum)
     %plot(log2(1:nsub),log2(v),'-*')
     plot(1:nsub,v,'k-')
     hold
     for n = 1:nsub
       plot([n n],[v(n)-seuil*stdasympSRD(n) v(n)+seuil*stdasympSRD(n)],'b-') ;   %SRD to compare
     end
     plot(1:nsub,mean(v)*line,'--')    %  add average of v's of the subseries as a dashed line
     
     %  Add CI's for entire trace on the left of variance plot
     stdasympSRD = sqrt(2/len)*vfull ;
     plot([.75 .75],[vfull-seuil*stdasympSRD vfull+seuil*stdasympSRD],'b-') ;
     %plot([.8,1:nsub],vfull+seuil*stdasympSRD.*CIline,'b-')  %  top line of CI
     %plot([.8,1:nsub],vfull-seuil*stdasympSRD.*CIline,'b-')  %  bottom line of CI
     plot(1:nsub,vfull*line,'k-')  %  add overall variance as a solid line
     
     V = axis;
     axis([ 0.5 nsub+0.5 V(3) V(4) ]); 
     hold off
     %xlabel('Subseries number')
     ylabel('Variances')
     %title('Variance of the subseries. The solid (dashed) horizontal line gives the overall (average) value.')
     grid
   end
end
if (lookH)
   if (printout)
     nplotnum = nplotnum+1;
     subplot(nplotnum)
  
     if (wantH)
       plot(1:nsub,H,'k-')
       HL = H-seuil*stdH;
       HR = H+seuil*stdH;
       fullHL = fullH-seuil*fullstdH;
       fullHR = fullH+seuil*fullstdH;
       hold on
       for n = 1:nsub
         plot([n n],[HL(n) HR(n)],'r-') ;   % add confidence intervals
       end
       plot([.75 .75],[fullHL fullHL],'r-') 
       plot(1:nsub,mean(H)*line,'--')    %  add average of H's of the subseries as a dashed line
     else
       plot(1:nsub,alpha,'k-')
       alphaL = alpha-seuil*stdalpha;
       alphaR = alpha+seuil*stdalpha;
       fullalphaL = fullalpha-seuil*fullstdalpha;
       fullalphaR = fullalpha+seuil*fullstdalpha; 
       hold on
       for n = 1:nsub
         plot([n n],[alphaL(n) alphaR(n)],'r-') ;   % add confidence intervals
       end
       plot([.75 .75],[fullalphaL fullalphaR],'r-')
       plot(1:nsub,mean(alpha)*line,'--')   
     end
     
     if (wantH)
       fprintf('For the full series, (j1,j2)= (%d,%d),  H = %4.2f,  cf = %7.4f,  Q= %5.4f \n\n',j1,j2,fullH,cf,Q)
       plot(1:nsub,fullH*line,'k-')
       %plot([.8,1:nsub],fullH+seuil*fullstdH.*CIline,'r-')  %  top line of CI
       %plot([.8,1:nsub],fullH-seuil*fullstdH.*CIline,'r-')  %  bottom line of CI
       V = axis;
       axis([ 0.5 nsub+0.5 min(HL)-0.2*abs(min(HL)-mean(H)) max(HR)+0.2*abs(min(HL)-mean(H)) ]);
       %xlabel('Subseries number')
       ylabel('H')
       title('Scaling parameter of the subseries. The solid (dashed) horizontal line gives the overall (average) value.')
     else
       fprintf('For the full series, (j1,j2)= (%d,%d),  alpha = %4.2f,   cf = %7.4f,  Q= %5.4f \n\n',j1,j2,fullalpha,fullcf,Q)
       plot(1:nsub,fullalpha*line,'k-')
       %plot([.8,1:nsub],fullalpha+seuil*fullstdalpha.*CIline,'r-')  %  top line of CI
       %plot([.8,1:nsub],fullalpha-seuil*fullstdalpha.*CIline,'r-')  %  bottom line of CI
       V = axis;
       axis([ 0.5 nsub+0.5 min(alphaL)-0.2*abs(min(alphaL)-mean(alpha)) max(alphaR)+0.2*abs(min(alphaL)-mean(alpha)) ]);
       ylabel('\alpha')
       title('Scaling parameter of the subseries. The solid (dashed) horizontal line gives the overall (average) value.')
     end
     grid
   end

   %%%%%%  Constancy of alpha test
   if nsub>1 
    siglevel = 0.05;
    sigma = stdalpha(1);  % all the std's are the same for same size blocks
    mean(alpha);
    %std(alpha)
    %std(alpha)^2
    %%%%% Calculate for general test
    V = ( (std(alpha))^2*(nsub-1)/nsub )*nsub/sigma^2 ;   % statistic = S^2 * m/sig^2 
    %% the crit level is passed back, it can be compared against siglevel in the calling program
    crit_level = 1 - chi2_CDF( V, nsub-1); % test probability corresponding exactly to the data
 
    %%%%% Calculate for two-level test form each value of m1, in [1,nsub-1]
    crit_level2 = ones(1,nsub-1);
    m1 = max(1,floor(cutpoints/4));   % trick for A4 section 4.3, only calculate the one needed : m/4
    %m1 = max(1,floor(cutpoints/2));   % trick for A4 section 4.2, only calculate the one needed : m/2
    % for m1 = 1:(nsub-1)
      xbar = mean(alpha(1:m1));
      ybar = mean(alpha((m1+1):nsub));
      V2(m1) = abs(xbar-ybar)/sigma;
      scaledsigma(m1) = sqrt( nsub/m1/(nsub-m1) );
      crit_level2(m1) = 2*gauss_CDF( -V2(m1), 0, scaledsigma(m1) );   % gauss_CDF can't take vector q
    %end
    %V2;

    if printout
      if wantH        % use this for legend placement
         textposition = max(HR);
      else
         textposition = max(alphaR); 
      end 
      C = chi2_CDFinv(1-siglevel,nsub-1);           % boundary of critical region 
      if  V < C
        fprintf('\n***** General constancy test for alpha cannot be rejected at significance level %4.3f \n', siglevel)
        Handle=text(nsub*0.8-1+0.1,textposition,['\bf Not Rejected, (',num2str(100*(crit_level),3), '% sig)']) ;
        %Handle=text(nsub*0.8-1+0.1,textposition,'\bf Not Rejected') ;
        set(Handle,'fontsize',12);
      else
        fprintf('\n***** General constancy test for alpha was rejected at significance level %4.3f \n',siglevel)
        Handle=text(nsub*0.8-1+0.1,textposition,['\bf Rejected, (',num2str(100*(crit_level),3), '% sig)'] ) ;
        %       Handle=text(nsub*0.8-1+0.1,textposition,'\bf Rejected') ;
        set(Handle,'fontsize',12); 
      end
      fprintf('      The test statistic was %4.2f corresponding to a critical level of %4.3f and probability of %4.3f \n',V,crit_level,1-crit_level)
      fprintf('      The critical region was to the right of %4.3f  \n', C)
    
      %% two-level output
      C2 = gauss_CDFinv(1-siglevel/2, 0, scaledsigma );    % critical regions get bigger with m1
      fprintf('\n      The two-level tests, beginning at m1 = 1,  gave critical levels of: \n')
      crit_level2
      fprintf('      Corresponding rejections ( 1=reject ) were:\n')
      rejections = (V2 >  C2);
      rejections

      hold off
    end  

   end  % constancy test

   %ans=input('Should I launch LDestimate on the trace? (hit return for yes)  ');
   %if  isempty(ans)
      %LDestimate(x,regu,1,j2,1)
   %end
end


if (lookcf)
   if (printout)
     nplotnum = nplotnum+1;
     subplot(nplotnum)
     plot(1:nsub,cf,'k-')
     sig_level = 5 ;  %set confidence level
     m = log( (cf.^4)./((stdcf.^2)+cf.*cf) )/2;
     v = log( (stdcf.^2)./cf./cf + 1 );
     cfL = exp(gauss_CDFinv(  sig_level/2/100,m,sqrt(v)));
     cfR = exp(gauss_CDFinv(1-sig_level/2/100,m,sqrt(v)));

     m = log( (fullcf^4)/((fullstdcf^2)+fullcf*fullcf) )/2;
     v = log( (fullstdcf^2)/fullcf/fullcf + 1 );
     fullcfL = exp(gauss_CDFinv(  sig_level/2/100,m,sqrt(v)));
     fullcfR = exp(gauss_CDFinv(1-sig_level/2/100,m,sqrt(v)));

     hold on
     for n = 1:nsub
       plot([n n],[cfL(n) cfR(n) ],'r-') ;   % add confidence intervals
     end
     plot([.75 .75],[fullcfL fullcfR],'r-')
     plot(1:nsub,mean(cf)*line,'--')    %  add average of cf's of the subseries as a dashed line

     plot(1:nsub,fullcf*line,'k-')
     %plot([.8,1:nsub],fullcf+seuil*fullstdcf.*CIline,'r-')  %  top line of CI
     %plot([.8,1:nsub],fullcf-seuil*fullstdcf.*CIline,'r-')  %  bottom line of CI
     V = axis;
     axis([ 0.5 nsub+0.5 min(cfL)-1.6*abs(min(cfL)-mean(cf)), max(cfR)+0.8*abs(min(cfR)-mean(cf)) ]);
     %xlabel('Subseries number')
     ylabel('c_f')
     title('Second scaling parameter of the subseries. The solid (dashed) horizontal line gives the overall (average) value.')
     hold off
     grid
     fprintf('For the full series, (j1,j2)= (%d,%d),  alpha = %4.2f,  cf = %7.4f,   Q= %5.4f \n\n',j1,j2,fullalpha, fullcf,Q)
   end
end

if (printout)
  xlabel('Block number')
  fprintf('\n*******************************************************************************************\n  \n')
end
