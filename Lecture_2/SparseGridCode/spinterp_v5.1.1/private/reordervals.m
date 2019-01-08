function z = reordervals(z, seq, fromindex, toindex)
% REORDERVALS   Reorder hierarchical surpluses of Chebyshev grid   
%    Private function that reorders the entries in z such 
%    that the dimensions which have a higher number of support nodes 
%    come first. (Internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Date   : January 24, 2006
% Version: 1.1

% Change log:
% V1.0   : July 7, 2004
%          Initial version
% V1.1   : January 24, 2006
%          Changed data types to operate on uint arrays

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

% This reordering results in a slight performance gain, since the
% inner loops always have a higher number of iterations compared to
% the outer loops.
	
if nargin < 3
	% Do reordering for full index array
	
	d = uint16(size(seq,2));
	if d > 1
		nlevels = uint32(size(seq,1));
		n = zeros(1,d,'uint32');
		index = uint32(1);
		for kl = 1:nlevels
			[newseq, order] = sort(seq(kl,:),2,'descend');
			npoints = uint32(1);
			for k = 1:d
				lval = seq(kl,k);
				level(k) = lval;
				if lval == 0
					n(k) = 1;
				else
					if lval < 3
						n(k) = 2;
						npoints = npoints * 2;
					else
						n(k) = 2^uint32(lval-1);
						npoints = npoints * 2^uint32(lval-1);
					end
				end
			end
			ztemp = z(index:index+npoints-1);
			ztemp = reshape(ztemp, n);
			ztemp = permute(ztemp, order);
			z(index:index+npoints-1) = ztemp(:);
			index = index + npoints;
		end
	end
else
	% Do the sparse indices version
	
	% Get the indices data from the structure (required here such that
	% the processing in the loops gets JIT-compiled).
	subGridPoints = seq.subGridPoints;
	subGridAddr = seq.subGridAddr;
	indicesNDims = seq.indicesNDims;
	indicesDims = seq.indicesDims;
	indicesAddr = seq.indicesAddr;
	indicesLevs = seq.indicesLevs;
	
	d = size(seq.forwardNeighbors,2);
	n = zeros(1,d,'uint32');

	currentindex = fromindex;
	nkl = uint32(1);
	
	index = uint32(1);
	while currentindex <= toindex
		% Get the number of grid points of the currently processed index
		npoints = subGridPoints(currentindex);
		ndims = indicesNDims(currentindex);
		addr = indicesAddr(currentindex);
	
		if ndims > 1
			[dummy, order] = ...
					sort(indicesLevs(addr:addr+uint32(ndims)-1),1,'descend');
			did = uint8(1);
			while did <= ndims
				lval = indicesLevs(addr);
				if lval < 3
					n(did) = 2;
				else
					n(did) = 2^uint32(lval-1);
				end
				addr = addr + 1;
				did = did + 1;
			end
			ztemp = z(index:index+npoints-1);
			ztemp = reshape(ztemp, n(1:ndims));
			ztemp = permute(ztemp, order);
			z(index:index+npoints-1) = ztemp(:);
		end
		index = index + npoints;
		currentindex = currentindex + 1;
	end
end
