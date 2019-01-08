function fullseq = spseq2full(seq)
% SPSEQ2FULL   Convert sparse sequence of index sets to full one
% (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : March 22, 2005

% Change log:
% V1.0   : March 22, 2005
%          Initial revision
% V1.1   : January 24, 2006
%          Minor update to generate uint8 instead of double array
           
% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

d = uint16(size(seq.forwardNeighbors, 2));
nlevels = uint32(size(seq.indicesNDims, 1));

fullseq = zeros(nlevels, d, 'uint8');
addr = uint32(1);
for k = 1:nlevels
	ndims = seq.indicesNDims(k);
	for did = 1:ndims
		fullseq(k,seq.indicesDims(addr)) = seq.indicesLevs(addr);
		addr = addr + 1;
	end
end
