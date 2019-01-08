function x = spgridgp(levelseq)
% SPGRIDGP  Compute grid points, Gauss-Patterson grid
%    X = SPGRIDGP(LEVELSEQ)  Computes the sparse grid points for
%    the given sequence of index LEVELSEQ. The coordinate value of
%    dimension i is stored in column i of the matrix X. One row of
%    matrix X represents one grid point.
%    (Internal function)
	
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
	
% Get the number of levels
nlevels = uint32(size(levelseq,1));

% Get the dimension
d = uint16(size(levelseq,2));

npoints = zeros(nlevels, 1, 'uint32');
dim = zeros(d,1, 'uint16');

% Compute number of points
totalpoints = uint32(0);
for k = 1:nlevels;
	ntemp = uint32(1);
	for l = 1:d
		ntemp = ntemp * 2^uint32(levelseq(k,l));
	end
	npoints(k) = ntemp;
	totalpoints = totalpoints + ntemp;
end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
x = 0.5*ones(totalpoints,d);
	
for kl = 1:nlevels
	c = {};
	ndims = uint16(0);
	for k = 1:d
		% compute the points, scaled to [0,1]
		lev = double(levelseq(kl, k));
		if lev == 0
			continue;
		else
			ndims = ndims + 1;
			ctemp = gpabsc(lev);
			c{ndims} = ctemp(1:2:end)';
		end
		dim(ndims) = k;
	end
	if ndims > 1
		[c{:}] = ndgrid(c{:});
	end
	for k = 1:ndims
		x(index:index+npoints(kl)-1,dim(k)) = c{k}(:);
	end
	index = index + npoints(kl);
end
