%---------------------------
%
% wtspec.m
% 
% PA DV 97-10-30
%
%
% rlistcoefdaub.m
%-------------------------------------

function [muj,nbj]=wtspec(appro,N,nbvoies) ;

nj = length(appro) ;
h1 = rlistcoefdaub(N) ;
nl = length(h1) ;
g1 = (-1).^(0:-1+nl).*fliplr(h1) ;
gg1 = fliplr(g1) ;
hh1 = fliplr(h1) ;

for j=1:nbvoies,
         convolue=conv(appro,gg1) ;
         decime=convolue(nl:2:nj) ;         %  decime always becomes empty near the end, 
         if length(decime) == 0
            break
         end
         muj(j) = mean(decime.^2) ;           %  generates a error here
         nbj(j)=length(decime) ;
         clear convolue decime
% ---compute the appro
         convolue  =conv(appro,hh1) ;
         appro = convolue(nl:2:nj) ;
         nj = length(appro) ;
         clear convolue 
end
%index=find(nbj> 2* N); %arbitraire
index = find(nbj >= 2 ) ;
muj=muj(index) ;
nbj=nbj(index) ;









