function ipder = dctdiffcheb(n,z,y)
% DCTDIFFCHEB  Compute the derivative of a univariate 
%    polynomial given by hierarchical surpluses at the
%    Cehbyshev-Gauss-Lobatto nodes) using FFT/IFFT. 
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : June 25, 2006

% Change log:
% V1.0   : Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

ninterp = size(z,2);
ipder = zeros(length(y),1);

% Extend grid values
if n == 3
	zext = [z(2,:); zeros(1,ninterp); z(1,:); zeros(1,ninterp)];
else
	zext = zeros((n-1)*2,ninterp);
	zext(2:2:end,:) = [flipud(z); z];
end

% Perform DCT via FFT
chebc = real(fft(zext))/(double(n)-1);

chebc(n,:) = chebc(n,:)/2;
chebc(1,:) = chebc(1,:)/2;

% Differentiate (in place)
for k = 1:ninterp
	c = chebc(n-1,k);
	chebc(n-1,k) = 2*double(n-1)*chebc(n,k);
	chebc(n,k)   = 0;
	l = double(n) - 2;
	while l >= 2
	  cprev = c;
		c = chebc(l,k);
		chebc(l,k) = chebc(l+2,k) + 2*l*cprev;
		l = l - 1;
	end
	chebc(1,k) = chebc(3,k)/2 + c;
end

% Normalize to interval [0,1]
chebc = chebc * 2;

% Evaluate using the Clenshaw recurrence formula
if ninterp ~= length(y)
  % Univariate case is slightly different; there is only
	% a single column of surplus values in z
	ninterp = length(y);
	for k = 1:ninterp
		c  = 0;
		c2 = 0;
		l  = n - 2;
		y1 = 2 * y(k) - 1;
		y2 = 2*y1;
		tmp = 0;
		while l >= 1
			tmp = c;
			c = y2*c - c2 + chebc(l+1);
			c2 = tmp;
			l = l - 1;
		end
		ipder(k) = y1*c - c2 + chebc(1);
	end
else
  % Process multiple columns from multivariate case
	for k = 1:ninterp
		c  = 0;
		c2 = 0;
		l  = n - 2;
		y1 = 2 * y(k) - 1;
		y2 = 2*y1;
		tmp = 0;
		while l >= 1
			tmp = c;
			c = y2*c - c2 + chebc(l+1,k);
			c2 = tmp;
			l = l - 1;
		end
		ipder(k) = y1*c - c2 + chebc(1,k);
	end
end
