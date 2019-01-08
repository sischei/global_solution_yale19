function ip = spinterpgp(d,z,y,levelseq,purgedata)
% SPINTERPGP  Variable degree polynomial interpolation (Gauss-Patt.)
%    IP = SPINTERPGP(N,D,Z,Y)  Computes the interpolated
%    values at [Y1, ..., YN] over the sparse grid at level N. Y
%    may be a double array for vectorized processing, with each row
%    representing one point. The sparse grid data must be given  as
%    an array Z containing the hierarchical surpluses (computed with
%    SPVALS). Note that the sparse grid is always normalized to
%    the unit cube [0,1]^D, i.e., if the weights have been computed
%    for a different domain, the values Y have to be rescaled
%    accordingly. (Internal function)
	
% Author : Andreas Klimke
% Version: 1.0
% Date   : December 22, 2007
	
% Change log:
% V1.0   : December 22, 2007
%          Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

ninterp = uint32(size(y,1));
ip = zeros(ninterp,1);

% Get the number of levels
nlevels = uint32(size(levelseq,1));

if ~isempty(purgedata), purge = true; else purge = false; end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
n = ones(d,1,'uint32');
order = ones(d,1,'uint16');

for kl = 1:nlevels
	npoints = uint32(1);
			
	lval = uint8(0);
	ndims = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		if lval > 0
			ndims = ndims + 1;
			n(ndims) = 2^uint32(lval+1)-1;
			npoints = npoints * 2^uint32(lval);
			order(ndims) = k;
		end
	end
	
	% Skip subgrids with all surpluses below droptol.
	if purge
		if purgedata(kl) == 0
			index = index + npoints;
			continue;
		end
	end
	
	if npoints > 1
		vals = z(index:index+npoints-1);
		nord = n(1:ndims);
		ip = ip + barypdstepgp(vals, nord, order(1:ndims), getgpnodes(nord), y, ...
				getgpbaryw(nord));
	else 
		% no interpolation necessary; we have just a single constant function
		ip = ip + z(index);
	end
	index = index + npoints;
end

