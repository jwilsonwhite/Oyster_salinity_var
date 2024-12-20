function [kmat] = kernmatSimp(x,Sal,Params,M,Type)

% Set up the integration mesh kernel using Simpsons
% For Oyster IPM

%adapted from Will from Easterling. Evenly spaced grid, now add weights so
%use Simpson's rule using Marissa's code.

y = x;
%this creates a vector (y) that is equal to x

[x,y] = meshgrid(x,y); % Matlab built in function
%x is an array (original x by original x) with each row equal to original
%vector x
%y is an array (original y by original y) with each column equal to
%original vector y
%X corresponds to size at time t
%Y corresponds to size at time t+1

% Make the kernel from the grid
kmat = mkkern(x,y,Sal,Params,M,Type);

kmat = max(0,kmat); 


