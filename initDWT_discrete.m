%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   initDWT_discrete.m                                  %
%                                                       %
%        D. Veitch   P.Abry                             %
%                                                       %
%   LYON 98-09-10                                       %
%   DV Melbourne  15/1/99                               %
%                    7/99                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function performs the initialization step of the DWT for a special kind
%    of process linked to the discrete time series X(k) supplied. 
%    In effect it enables the DWT to be used to study such intrinsically discrete time series.
%    From the X(k) supplied a band limited cts time process ~X(t) is defined, with X(k) = ~X(k)
%    via the sinc expansion. This function calculates the approximation coefficients for ~X(t):
%    
%           a(k) = integral^infty_infty  ~X(t) phi(t-k) dt   
%
%    where phi(t) is the scale function of the MRA underlying the DWT. The integrals
%           I(m) = integral^infty_infty  sinc(t+m) phi(t) dt  have been precalculated in
%    calc_initfilterDWT_discrete.m  and so this routine can calculate the a(k) simply via
%    the discrete convolution:
%                                 a(k) = (X*I)(k)   
%
%--- Routines called Directly :
%
% calc_initfilterDWT_discrete.m
%      " function [I,nI] = calc_initfilterDWT_discrete(regu,seuil) "
%
% initfilterDaub_DWT_discret.m
%      " function [I,nI] = initfilterDaub_DWT_discrete(regu) "
%
%  Input:   data:  the discrete time series  X(k),  k in [1,2,... n] 
%           regu:  (regularity) =  number of vanishing moments (of the Daubechies wavelet).
%           lengthIwant:   the desired filter length.  If =0, the automated length selection method
%                       is chosen:  symmetric truncation if necessary to 1/8th of the data length. Otherwise
%                       the requested length is used if possible (can't exceed the stored length).
%                       It is not allowed however to exceed the data length.  Warning messages are
%                       printed in these cases.
%           printout:  if = 1 , have plots and information on the filter and the approximation series
%                         = 0 , no output except the warning messages if applicable.
%
%  Output:  appro: The approximation coefficents at level zero (set to correspond to deltat=1)
%           kfirst: The time t=k which the first element of appro corresponds to.     
%           klast:  The time t=k which the last  element of appro corresponds to.
%         For both of these time is on the scale of t in ~X(t), which is that of k in X(k).
%

function [appro,kfirst,klast] = initDWT_discrete(data,regu,lengthIwant,printout) ;

n = length(data);

%--- Compute or load the sequence I(m) involved in the initialisation
  %--- If not DaubechiesN with N=1 or 2 or .... 10, then need to calculate the filter
  %fprintf('Beginning calculation of the I(m)....  ')
  %[I,nI] = gen_initfilterDWT_discrete(regu,precision) ;
  %fprintf('completed, length of filter is %d \n\n',length(I))
  
  %--- Load the precalculated filter if DaubechiesN with N=1 or 2 or .... 10.
  [I,nI] = initfilterDaub_DWT_discrete(regu);  
  lenI = length(I);


%--- Determine the length of the filter.
%    The length can be specified directly by lengthIwant,  or chosen automatically
%    whereby if data is too short, it is truncated to avoid generating an approximation series
%    which is too short. 
%    It is not allowed to be longer than the data, and if a filter length is requested which is
%    longer than the stored filter, then the whole filter is used (same as if lengthIwant=0).

  if lengthIwant==0  %  code for the automatic choice
    % ***  Truncate if necessary so the filter is no more than one eighth as long as the data
    % Note that this approach gives a variable precision dependent on N (regu), however when
    % dealing with short series, controlling the Length is paramount. The precision is indicated 
    % by printing the end values of the filter.  
    if  lenI > n/8
      if printout 
        fprintf('Data is small, symmetrically truncating the filter to 1/8th of data length.\n');
      end
      radius = floor(n/16);
      % lengthIwant = 19;                 % if want to fix the length manually, do it via this here
      % radius = ( lengthIwant - 1 )/2;
      leftside  = max(1,   nI-radius);
      rightside = min(lenI,nI+radius);
      I = I(leftside:rightside);          % will normally be symmetric:  length 2*radius+1
      lenI = length(I);
      nI = radius + 1;
    end
  else   %  use the input value, with a check 
    if  lengthIwant > n 
      fprintf('Requested filter length exceeds data length of: %d, truncating to data length\n',n)
      lengthIwant = min(lengthIwant,n) ;
    end
    if  lengthIwant > lenI      %  just use the stored filter, ignore lengthIwant
      fprintf('Truncating requested filter length to length of stored filter: %d \n',lenI)
    else  % truncate the stored values symmetrically to the requested length
      radius = floor(( lengthIwant - 1 )/2);   % ensure an integer radius
      leftside  = max(1,   nI-radius);     
      rightside = min(lenI,nI+radius);
      I = I(leftside:rightside);           % will normally be symmetric:  length 2*radius+1
      lenI = length(I);                    % can be asymmetric is symmetric truncation finds an end. 
      nI = radius + 1;
    end
  end

 
  % plot the filter, plus two closeups
  if printout
    fprintf('Filter used is  %d elements long.  The %dth element corresponds to m=0\n',lenI,nI);
    fprintf('The end elements are:  ( I(1), I(%d) ) = ( %f, %f ) \n',lenI,I(1),I(lenI) );
    figure(2)
    clf
    index = (1-nI):(-nI+lenI);   % calculate indices for plotting
    subplot(311) ; plot(index,I,'*')    ; grid
    title('non-negligible filter coefficients I(m)')
    width = min(3, floor(lenI/2));
    subplot(312) ; plot(index((nI-width):(nI+width)), I((nI-width):(nI+width)),'*') ; grid
    title('close up centred on m=0')
    subplot(313) ; plot(index(max(1,lenI-5):lenI), I(max(1,lenI-5):lenI),'*') ; ; grid
    title('Last 6 values')
  end

%--- perform the initialization of the continuous time process ~X associated to the data 
%   With the scale function having support on t in [0,lenphi] = [0,2N-1] for DaubechiesN,
%    and assuming that ~X(t) in known perfectly for t in [1,n], the allowed coefficients
%    are k = 1, 2, ... n-(lenphi)=n-2N+1  for DaubechiesN.
%   However the ~X(t) are not known perfectly, and sinc drops off quite slowly.
%    The decision of how many to exclude at the edges has already been made 
%    implicitly in the calculation of I(m).  The corresponding polluted coefficients of 
%    appro will be  trimmed here.  On the right typically more will be trimmed
%    due to this effect than the other (this is checked for).
%    (note that with k=1 (impossible value), I(0) is aligned with X(1) )
%                    k=n                     I(0)                 X(n)
appro = conv(data,I) ;

% trim polluted coefficients due to conv going over the edges
appro = appro(lenI:n) ;     % length is  n-lenI+1

% determine the k values the ends of the unpolluted series correspond to
lengthleft = nI-1;          % number of elements in I to the left of m=0
lengthright = lenI - nI;    %                                right
kfirst = 1 + lengthright;   % left and right are reversed since one convolves
klast  = n - lengthleft;    % length of range (kfirst,klast) is also  n-lenI+1  as reqd
if klast > n - 2*regu + 1
  if printout 
    fprintf('klast limited by support of scale function\n');
  end
  klast = n - 2*regu + 1;
  appro = appro(1: klast - kfirst +1);
end

% compare data and approximation series.
if printout
  fprintf('Data is from k=1 to k=%d.\n',n)
  fprintf('given phi edge effects, max possible range for appro is k=1 to k=%d.\n',n-2*regu+1)
  fprintf('with sinc edge effects, max possible range for appro is k=%d to k=%d, this is the output.\n',kfirst,klast)
  figure(1)
  clf
  plot(1:n,data,'b-');
  hold on
  plot(kfirst:klast,appro,'r-');
  title('Discrete data is in blue, the approximation coefficients in red')
  hold off
end
