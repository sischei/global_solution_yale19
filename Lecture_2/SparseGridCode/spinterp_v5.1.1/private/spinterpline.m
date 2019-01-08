function [ip,ipgrad] = spinterpline(z,p,xdir,x)
% SPINTERPLINE   Evaluate the sparse grid interpolant along a line
%    IP = SPINTERPLINE(Z,P,XDIR,X)  Evaluate the D-dimensional 
%    interpolant at a point Y = P + XDIR*X, where P and XDIR are 
%    Dx1 vectors, and X is a scalar. This function can also be
%    called in vectorized form to evaluate multiple points along
%    the line at once. In this case, X must be a row or column
%    vector.
%
%    [IP,IPGRAD] = spinterpline(...)  Returns also the directional
%    derivative IPGRAD at the point Y = P + XDIR*X with respect to
%    X.
%
%    Example:
%       f = inline('x.^2 + y.^2 - 2.*z');
%       z = spvals(f,3);
%       p = [0 0 0];
%       xdir = [1 1 1];
%       x = linspace(0,1,6);
%       f_interp = spinterpline(z, p, xdir, x)

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

% Make p a column vector
p = p(:);   

% Make xdir a row vector
xdir = xdir(:)';

% Check if size of x is valid
xsize = size(x);
if length(xsize) > 2 || ( xsize(1) > 1 && xsize(2) > 1 )
	error('Matlab:spinterp:badopt','X must be either a row or a column vector.');
end	

% Prepare points to evaluate
nx = length(x);
p = repmat(p,[1,nx]);
if xsize(1) > 1
	% x is a column vector
  x = p' + x * xdir(:)';
	xcell = num2cell(x,1);
else
	x = p + xdir(:) * x;
	xcell = num2cell(x,2);
end

% Evaluate and return results
if nargout > 1
	[ip,ipgradvec] = spinterp(z,xcell{:});
	ipgrad = zeros(xsize);
	for k = 1:nx
		ipgrad(k) = xdir * ipgradvec{k}; 
	end
else
	ip = spinterp(z,xcell{:});
end
