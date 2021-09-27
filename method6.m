%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   method6.m                                           %
%                                                       %              
%        D. Veitch   P.Abry                             %
%                                                       %
%   Melb 1/5/2000                                       %
%  DV: Melb 1/6/2002, method fine tune during write up, %
%                     support for fac passing           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Implements method 6 of selecting j1*, the lower cutoff for LRD, based on a statistic 
%   passed to it, given as a function of j1.
%
%   Method:   Define a non-decreasing zone beginning from j1.  
%             Within this zone find the LASt time that an improvement (ratio of Q's) larger
%             than 'fac' was found.  
%        Difference from method4:  then add one to it, as long one stays within the nondec zone
%             If none found, or never decreases,  j1opt=1 .
%
%   Input:    Qmat:  a matrix, each row is for a j2 fixed,  indexed by j1, the measure of Quality 
%                    of the estimate given j1
%             j2vec: vector (possibly corrected) of j2 values through from newchoosej1
%  
%   Output:  j1opt:  the j1* values chosen, one per row of the matrix (different j2 values), as a vector.
%             logQ:  a modified log of Q vector suitable for plotting  (lower bounded at -10)
%                     (-10 also appears in 'empty' positions in the matrix (since j2 varies..))
%       endofnodec:  the location of the end of the no decreasing search zone (one value per j2)
%----------------------------------------------------------------------------------------------
function [j1opt,logQ,endofnodec] = method6(Qmat,j2vec)
%function [j1opt,logQ,endofnodec] = method6(Qmat,j2vec,realfac)        % trick for article2, pass fac

%%%%% constants used in the method
lowerbound = 1.e-10;     % value is arbitrary, but shouldn't introduce dependence as so small
fac = 10;                % the 'improvement factor', arbitrary mult. factor, search for larger jumps than this.
%fac = realfac;           % trick for article2, pass fac

%%%%% Process input values
Qmat(Qmat<lowerbound) = lowerbound;    % lower bound the Q's to make it easier to plot, and avoid zero value.
logQ = log10(Qmat);                    % store log values for returning, and the heuristic phase...
lenj2 = size(Qmat);  lenj2 = lenj2(1); % determine number of rows (different j2 values)

%%%%% plot log Q values (diagnostic)
%figure(42)
%plot(logQ');grid

%%%%% Apply the heuristic method to the statistic Q (although for convenience work with logQ)
%%%--- loop over  j2  values 
for k = 1:lenj2 
   lQ = logQ(k,:);     %  log Q for this j2
   %lQ = [-1 -1.5  -1.4  -1.5 -1 -.05 -.04 -.05 ];           % testvalues   nasty case , decrease immediately
   %lQ = [-10 -10 -9 -4 -2 -1.4 -1.3 -1 -.05 -.04 -.03 ];    % testvalues, nasty case,  never decreases
   %lQ = [-10 -10 -9 -4 -2 -1.4 -1.5 -1 -.05 -.04 -.05 ];    % testvalues, typical case
   %figure(44); plot(lQ);grid

   %%%%% apply the algorithm
   factors = diff(lQ);                % this generates the factors, as in the log domain
   signchanges = find(factors<0);     % If want to change to INC instead of non-dec do it here
   if length(signchanges)==0
      endofnodec(k) = length(lQ);           % never decreases! 
   else
      endofnodec(k) = signchanges(1);       % find the end of the initial nondecreasing zone  
   end
   
   if endofnodec(k) == 1     % if have no choice, choose j1opt=1, even if bad. 
     j1opt(k) = 1;
   else                   % find Last time in endofnodec zone when improvement bigger than fac
     bigjumps = find(factors(1:endofnodec(k)-1)>= log10(fac));  %  which jumps in Q in range are bigger than fac
     if length(bigjumps) == 0                                %  no big jumps, so stay at j1=1
       j1opt(k) = 1;
     else                                                    %  find last one
       j1opt(k) = bigjumps(length(bigjumps))  + 1 ;          %  plus 1 because took diff
     end
     %j1opt(k) = min( j1opt(k)+1, endofnodec(k));            % refuse to enter dec zone 
     J1opt(k) = j1opt(k) + 1;                                % original 4a + 1, seems to be better 
     j1opt(k) = min(j1opt(k),j2vec(k)-2);                    % do not exceed maximum possible j1 value for this j2
  end
   %j1opt(k) = j1opt(k) + 1;            % add 1 in all cases, makes model 9 worse but more coherent
end   %%  loop over j2

















