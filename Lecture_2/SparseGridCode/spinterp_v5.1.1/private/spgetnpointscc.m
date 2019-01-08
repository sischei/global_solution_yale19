function [totalpoints, npoints] = spgetnpointscc(levelseq);
% SPGETNPOINTSCC   Compute the number of subgrid points (CC-Grid)
%    (internal function)
	
% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : February 1, 2006
	
% Change log:
% V1.0   : February 1, 2006
%          Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

nlevels = uint32(size(levelseq, 1));
d = uint16(size(levelseq, 2));
npoints = zeros(nlevels,1,'uint32');

totalpoints = uint32(0);
% Compute the number of points per subdomain of new levels
for kl = 1:nlevels
	np = uint32(1);
	lval = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		if lval == 0 % do nothing
		elseif lval < 3
			np = np * 2;
		else
			np = np * 2^uint32(lval-1);
		end
	end
	npoints(kl) = np;
	totalpoints = totalpoints + np;
end	
