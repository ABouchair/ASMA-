%ffteur.m
%
%   FT of signal x with output n points- the input is padded to n points if necessary, and the output
%   of a the FFT of a vector of length n is also of length n
%
%   The freq range normalisation is decided manually, as determined by fe
%   ex  fe = 1, n=8   nu = [ -1/2 -3/8 -1/4 -1/8 0 1/8 1/4 3.5/8 ] 
%
%    so frequency output is:  f = [ -fe/n .. 0 .. (n-1)fe/n ]

function [f1,f2]=ffteur(x,n,fe);

f1=fft(x,n);       % FFT of signal x with n points in the output, in normalised freq nu in [0,1/n)
                   % normalised frequencies  from 0 to (n-1)/n returned in elements 1 to  n
                   % DFT is symmetric as signal is real
f1=f1/fe;          % renormalize the frequencies to f = [0, fe/n)  
f1=fftshift(f1);   % convert to a form where f = 0 is the middle term, ie  [-fe/2, fe/2)
                   % (since n is even, this implies that nu=1/2 (f = fe/2) is dropped)  
f2=(-n/2:1:(n-1)/2)/n*fe;        %  the corresponding freq range
