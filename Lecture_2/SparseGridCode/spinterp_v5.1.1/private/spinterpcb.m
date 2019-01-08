function ip = spinterpcb(d,z,y,levelseq,purgedata)
% SPINTERPCB  Variable degree polynomial interpolation (Chebyshev)
%    IP = SPINTERPCB(N,D,Z,Y)  Computes the interpolated
%    values at [Y1, ..., YN] over the sparse grid at level N. Y
%    may be a double array for vectorized processing, with each row
%    representing one point. The sparse grid data must be given  as
%    an array Z containing the hierarchical surpluses (computed with
%    SPVALS). Note that the sparse grid is always normalized to
%    the unit cube [0,1]^D, i.e., if the weights have been computed
%    for a different domain, the values Y have to be rescaled
%    accordingly. (Internal function)
	
% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.2
% Date   : February 3, 2006
	
% Change log:
% V1.0   : July 5, 2004
%          Initial release.
% V1.1   : January 24, 2006
%          Changed data types to operate on uint arrays
% V1.2   : February 3, 2006
%          Added droptol processing.

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
level = ones(d,1,'uint8');

for kl = 1:nlevels
	npoints = uint32(1);
			
	lval = uint8(0);
	ndims = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		level(k) = lval;
		if lval == 0
			n(k) = 1;
		else
			ndims = ndims + 1;
			n(k) = 2^uint32(lval)+1;
			if lval < 3
				npoints = npoints * 2;
			else
				npoints = npoints * 2^uint32(lval-1);
			end
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
		[t, order] = sort(level,1,'descend');
		order = uint16(order(1:ndims));
		nord = n(order);
		vals = z(index:index+npoints-1);
		ip = ip + barypdstepcb(vals, nord, order, getchebnodes(nord), y);
	else 
		% no interpolation necessary; we have just a single constant function
		ip = ip + z(index);
	end
	index = index + npoints;
end
