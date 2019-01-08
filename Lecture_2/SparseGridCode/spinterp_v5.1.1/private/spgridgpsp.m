function x = spgridgpsp(seq, fromindex, toindex)
% SPGRIDGPSP  Compute grid points, Gauss-Patterson (sparse indices)
%    X = SPGRIDGPSP(SEQ, FROMINDEX, TOINDEX)  Computes the sparse
%    grid points for the given sequence of index sets SEQ, starting
%    with index FROMINED up to index TOINDEX. The coordinate value
%    of dimension i is stored in column i of the matrix X. One row
%    of matrix X represents one grid point. This routine works with
%    the sparse index set data generated with SPGETSEQSP.
%    (Internal function)
%	
% See also SPGETSEQSP.
	
% Author : Andreas Klimke
% Version: 1.0
% Date   : November 18, 2007

% Change log:
% V1.0   : November 18, 2007
%          Initial version

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

if nargin < 2
	if isfield(seq, 'currentindex')
		fromindex = seq.currentindex;
	else
		fromindex = uint32(1);
	end
	toindex = size(seq.indicesNDims,1);
end

% Get the number of levels
nlevels = toindex - fromindex + 1;

d = uint16(size(seq.forwardNeighbors,2));

% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);

totalpoints = uint32(sum(seq.subGridPoints(fromindex:toindex)));
x = 0.5*ones(totalpoints,d);
dim = zeros(d,1,'uint16');

currentindex = fromindex;
while currentindex <= toindex
	ndims    = seq.indicesNDims (currentindex);
	addr     = seq.indicesAddr  (currentindex);
	npoints  = seq.subGridPoints(currentindex);
	
	c = cell(double(ndims),1);
	k = uint8(1);
	while k <= ndims
		lev = double(seq.indicesLevs(addr));
		dim(k) = seq.indicesDims(addr);
		ctemp = gpabsc(lev);
		c{k} = ctemp(1:2:end)';
		addr = addr + 1;
		k = k + 1;
	end
	if ndims > 1
		[c{:}] = ndgrid(c{:});
	end
	k = uint8(1);
	while k <= ndims
		x(index:index+npoints-1,dim(k)) = c{k}(:);
		k = k + 1;
	end
	index = index + npoints;
	currentindex = currentindex + 1;
end
