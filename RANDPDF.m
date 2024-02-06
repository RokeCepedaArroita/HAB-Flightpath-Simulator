function x = RANDPDF(p,px,dim)
% RANDPDF
%   Random numbers from a user defined distribution
%
% SYNTAX:
%   x = randpdf(p, px, dim)
%       randpdf(p, px, dim)
% 
% INPUT:
%   p   - probability density,
%   px  - values for probability density,
%   dim - dimension for the output matrix.
%
% OUTPUT:
%   x   - random numbers. Run function without output for some plots.
%
% DESCRIPTION:
%   x = randpdf(p, px, dim) returns the matrix of random numbers from
%   probability density distribution defined in p and px. p are the density
%   (the y axis) and px are the value (the x axis) of the pdf. p and px
%   must be of the same length.
%   dim define the output matrix dimensions, for example dim=[100 3] define
%   the 100x3 two dimensional matrix with 300 random numbers.


% vectorization and normalization of the input pdf

px = px(:);
p = p(:)./trapz(px,p(:));

% interpolation of the input pdf for better integration

pxi = linspace(min(px),max(px),10000)';
pi = interp1(px, p, pxi,'linear');

% cumulative distribution function for input PDF

cdfp = cumtrapz(pxi,pi);

% finding the parts of cdf parallel to the X axis 

ind = [true; not(diff(cdfp) == 0)];

% and cut out the parts

cdfp = cdfp(ind);
pi = pi(ind);
pxi = pxi(ind);

% generating the uniform distributed random numbers

uniformDistNum = rand(dim);

% and distributing the numbers using cdf from input pdf

userDistNum = interp1(cdfp,pxi,uniformDistNum(:)','linear');

% making graphs if no output exists

x = reshape(userDistNum,dim);

end