function w = spquadwccsp(seq,w1d,id)
% SPQUADWCCSP  Compute quadrature weights, Clenshaw-Curtis grid 
%    W = SPQUADWCCSP(LEVELSEQ,W1D,ID)  Computes the quadrature 
%    weights for the given sequence of index sets LEVELSEQ. One 
%    row of column vectow W represents one weight. 
%    W1D are the integration weights in 1D, ID are the indices
%    where the weights of an according level start. 
%    This routine does essentially the same as SPQUADWCC, but 
%    operates over sparse index sets for increased efficiency, 
%    especially in case of higher dimensions D.
%    (Internal function)
%
% See also SPQUADWCC.

% Author : Andreas Klimke
% Version: 1.0
% Date   : September 13, 2007

% Change log:
% V1.0   : September 13, 2007
%          Initial version

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

indicesNDims  = seq.indicesNDims;
indicesLevs   = seq.indicesLevs;
subGridPoints = seq.subGridPoints;

% Get the number of levels
nlevels = uint32(length(indicesNDims));
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
% currentindex is the index of the currently processed multiindex
% entry. 
currentindex = uint32(1);
addr = uint32(1);

% Init weights vector
w = zeros(sum(subGridPoints),1);

while currentindex <= nlevels
	lval = 0;
	wval = 1;
	ndims = indicesNDims(currentindex);
	npoints = subGridPoints(currentindex);
	
	did = uint8(1);
	while did <= ndims
		lval = indicesLevs(addr);
		if lval < 3
			wval = wval * 1/4;
		else
			wval = wval * 1/(2^double(lval));
		end
		addr = addr + 1;
		did = did + 1;
	end

	w(index:index+npoints-1) = wval;
	
	index = index + npoints;
	currentindex = currentindex + 1;
end
