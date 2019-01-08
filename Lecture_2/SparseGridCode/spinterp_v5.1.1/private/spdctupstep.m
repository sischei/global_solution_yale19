function ip = spdctupstep(z, oldlev, newlev, dims)
% SPDCTUPSTEP  Step of polynomial interpolation at the Chebyshev-
%    Gauss-Lobatto nodes at grid nodes using DCTT 
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : February 8, 2006

% Change log:
% V1.0   : Initial release.
% V1.1   : February 8, 2006

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

d  = uint8(length(oldlev));
d2 = uint8(length(newlev));
n  = zeros(1,d,'uint32');

% Compute number of surpluses and number of chebyshev nodes in each dimension
for k = 1:d
	lval = oldlev(k);
	if lval < 3
		n(k) = 2;
	else
		n(k) = 2^uint32(lval-1);
	end
end

for k = 1:d
  % Check if there is anything to do in this dimension
  lval = oldlev(k);
	if newlev(k) > lval
		nact = n(k);

		% Make z-array d-dimensional
		if d > 1
			z = reshape(z, n);
			switchdims = [];
			if k > 1
				% Move current interpolation dimension to first dimension
				switchdims = [k 1:k-1 k+1:d];
				z = permute(z, double(switchdims));
				n = n(switchdims);
				dims(1:d) = dims(switchdims);
			end
		end
		
		% Reshape to 2D-array
		if d > 1
			ncols = prod(double(n(2:d)));
		else
			ncols = 1;
		end
		ztmp = reshape(z, [nact ncols]);
		% Pad with zeros
		if lval == 1
			z = zeros(3,ncols);
			z([1 3],:) = ztmp;
		else
			z = zeros(2^uint32(lval)+1,ncols);
			z(2:2:end,:) = ztmp;
		end
		% Update the value of n
		n(1) = 2^uint32(newlev(k)-1);
		% Perform DCT
		z = dctupsample(z, double(2^uint32(newlev(k))+1));
	end
end

% Now process the new dimensions. Here, we can just replicate the existing
% array.
nnew = [];
if d2 > d
	nnew = zeros(1,d2-d);
	for k = 1:d2-d
		lval = newlev(k+d);
		if lval < 3
			nnew(k) = 2;
		else
			nnew(k) = 2^uint32(lval-1);
		end
	end
	z = repmat(z(:), [1 nnew]);
end

% Finally, reshape z to its final form
if d > 1
	% Make array z d-dimensional
	z = reshape(z, [n nnew]);
end
if d2 > 1
  % Permute z to be sorted by dimension
  if d2 == 2
		if dims(1) > dims(2)
			z = permute(z, [2 1]);
		end
	else
		[dims, id] = sort(dims);
		dims = 1:d2;
		dims = dims(id);
		z = permute(z, double(dims));
	end
	z = z(:);
end
ip = z;
