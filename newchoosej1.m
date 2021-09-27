%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   newchoosej1.m                                       %
%                                                       %              
%        D. Veitch   P.Abry                             %
%                                                       %
% New version with separate statistic/heuristic phases, %
% DV 2/5/2000                                           %
% DV, Melb 18 May 2002  Add cheat for article 2         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    This function calculates a goodness of fit statistic for each choice of j1, then 
%    chooses an optimal j1 based on the value where sharp improvements are observed to
%    enter a stable region, if any.
%
%    Different heuristic methods may are implemented, using different 'method#.m' routines.  
%    For each a vector is returned corresponding to the vector of j2 values passed.
%
%    This routine need to be passed the muj, which can be conveniently accessed in LDestimate 
%    with the calling parameter  calcj1 = 1  .
%
%  Usage:  [j1opt,Qmat] = newchoosej1(regu, nj, muj, printout, 1:length(nj)))
%
%--- Routines called Directly :
%
%   regrescomp.m
%      " function [alphaest,cfCest,cfest,Cest,Q,Valpha,VcfC,CoValphacfC,Vcf,CoValphacf,unsafe]  
%                                                =  regrescomp(regu, nj, muj,j1,j2,printout) "
%   method#.m
%      " function [j1opt,logQ,endofnodec] = method#(Qmat,j2vec) "   # = 4 or 6
%
%
%  Input:   regu:  (regularity) =  number of vanishing moments (of the Daubechies wavelet).
%           nj:    the number of coefficients at scale j.
%           muj:   the vector[1..scalemax] of the average of the squared wavelet coefficients
%
%           printout: 1:  - then the graph of log_10( Q ) is plotted against j1 on FIGURE 5
%                         - the values of Q and the chosen j1 are printed 
%                     0:  - nothing is printed or plotted.
%
%           j2vec:  an arbitrary vector of j2 values
%
%
%  Output:  j1opt:     the chosen value of j1 according to the heuristic method chosen 
%           Qmat:      the matrix of Q values,  j1 = 1,2, ... j2  for each j2 in the input vector
%
%  the names of the chosen heuristics appear here:
%-----------------------------------------------------------------------------------------
%function [j1opt, Qmat, j1m5] = newchoosej1(regu, nj, muj, printout, j2vec) % cheat for article 2 
function [j1opt, Qmat] = newchoosej1(regu, nj, muj, printout, j2vec)

%%  Fix up possible bad input j2 values
maxj = length(nj);
if (min(j2vec)<3 | max(j2vec)>maxj)  
   j2vec = maxj;    % ensure j2 values are in the right range
end
lenj2 = length(j2vec);
maxj2 = max(j2vec);

%%  Internal parameters
fsize = 14;      % set font size (20 for papers)

%%%%% initialize
Qmat = zeros(lenj2,(maxj-2));         % each row of the Q matrix is for a j2 fixed.

%% selection of method
method = 6;


%%%%% PHASE 1:  calculate the quality statistic, here the original 'Q' probability values
for k = 1:lenj2      % loop over  j2  values
   j2 = j2vec(k);
   for j1= 1:j2-2    % loop over  j1 values 
     [alphaest,cfCest,cfest,Cest,Q] =  regrescomp(regu, nj, muj,j1,j2,0);   % calculate Q
     Qmat(k,j1) = Q ;  % store Q value
   end
end


%%%%% PHASE 2:  apply the heuristic(s) to the statistics, to choose the optimal value
if method==4
   [j1opt,logQ,endofnodec] = method4(Qmat);
elseif method==6
   [j1opt,logQ,endofnodec] = method6(Qmat,j2vec);
elseif method==5
   [j1opt,logQ,endofnodec] = method5(Qmat,j2vec);
end

%j1m5 = method5(Qmat,j2vec);    % calculate method 5 result also, passed back for article 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qmat = logQ;   % cheat for article2


if (printout)
  fprintf('\n******************  In "newchoosej1"  *****************************************************\n  \n');
  fprintf(' j2   j1*  j1= 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16\n\n')
  
  for i = lenj2:-1:1
    j2 = j2vec(i);
    fprintf('%3d:  %3d      ',j2,j1opt(i));
    for j1= 1:j2-2 
      fprintf('%3.2f ',Qmat(i,j1))
    end  
    fprintf('\n');
  end
  
  figure(10);
  hold off
  for i = lenj2:-1:1
    j2 = j2vec(i);
    j1= 1:j2-2;
    plot( j1, logQ(i,j1), 'r--', 'LineWidth',2 );         %  plot the Q values
    hold on
    %plot(j1(1:min(j2-2,endofnodec)),logQ(i,1:min(j2-2,endofnodec)), 'k-', 'LineWidth',2);  % nasty fix up, wait SERC files 
    plot(j1(1:endofnodec),logQ(i,1:endofnodec), 'k-', 'LineWidth',2 );  % plot first increasing part
    plot( j1opt(i),logQ(i,j1opt(i)), 'kd','markersize',10,'MarkerFaceColor','k') 
  end
  
  set(gca,'FontSize',fsize-2)
  xlabel('j_1');
  ylabel('log_{10} [ Q( j_1 ) ]');
  if length(j2vec)>1
     title(['Goodness of fit  Q( j_1) (for different j_2),  N= ',num2str(regu), ...
     '   (symbol gives j_1^* using method ',num2str(method),' )']);
  else
     title(['Goodness of fit  Q( j_1),  N= ',num2str(regu), ... 
     '   (symbol gives j_1^* using method ',num2str(method),' )']);
  end
  
  
  grid;
  %  add lines to remind us of the size of things... (though we try to avoid a heuristic that depends
  %  on the magnitude of Q itself
  plot( 0:maxj2-1 , -1*ones(1,maxj2), 'r-')            %  horizontal line at Q = 0.1
  plot( 0:maxj2-1 , log10(0.05)*ones(1,maxj2), 'r-')   %  horizontal line at Q = 0.05
  plot( 0:maxj2-1 , -2*ones(1,maxj2), 'r-')            %  horizontal line at Q = 0.01
  
  %% modify details of plot
  hhh=get(gca,'xlabel');
  set(hhh,'fontsize',fsize);
  hhh=get(gca,'ylabel');
  set(hhh,'fontsize',fsize);
  V = axis;
  axis([ 0 maxj2-1  min(-2,V(3)) 0 ])
  
  hold off
  
  fprintf('\n*******************************************************************************************\n  \n')
end












