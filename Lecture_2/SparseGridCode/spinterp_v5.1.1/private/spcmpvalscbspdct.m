function ip = spcmpvalscbspdct(d,z,y,seq,fromindex,toindex)
% SPCMPVALSCBSPDCT   Compute surpluses using DCT. Otherwise,
%    does the same as SPCMPVALSCBSP
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : February 14, 2006
	
% Change log:
% V1.0   : February 14, 2006
%          Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

% Allocate memory for the resulting surpluses
ninterp = size(y,1);
ip = zeros(ninterp,1);

% Get the first new active index
currentindex = fromindex;

% Compute total number of new indices that require computation of
% new surpluses
nnewlevels = toindex - currentindex + 1;

% Get the indices data from the structure (required here such that
% the processing in the loops gets JIT-compiled).
subGridPoints = seq.subGridPoints;
subGridAddr = seq.subGridAddr;
indicesNDims = seq.indicesNDims;
indicesDims = seq.indicesDims;
indicesAddr = seq.indicesAddr;
indicesLevs = seq.indicesLevs;

backward = seq.backwardNeighbors;
if isfield(seq, 'forwardAddr')
	forwardAddr = seq.forwardAddr;
	directForward = false;
else
	directForward = true;
end
forward = seq.forwardNeighbors;
	
% Allocate temporary arrays
level = zeros(1,d,'uint8');
newlevel = zeros(1,d,'uint8');
oldlevel = zeros(1,d,'uint8');
dimvec = ones(1,d,'uint16');
backvec = ones(d,1,'uint32');
olddimvec = ones(1,d,'uint16');
newdimvec = ones(1,d,'uint16');
newlevel2 = ones(1,d,'uint8');
down = ones(d,1,'uint8');

% kstart indicates the current start address of a subgrid of points
% of the new index sets.
kstart = uint32(1);

nkl = uint32(1);
while nkl <= nnewlevels
	% Get the number of grid points of the currently processed index
	npoints = subGridPoints(currentindex);
	
	% The number of subgrids contributing to the hierarchical
  % surpluses of current index.
	nsubgrids = uint32(1);
	
	newaddr = indicesAddr(currentindex);
	ndims = indicesNDims(currentindex);
	did = uint16(1);
	while did <= ndims
		dim = indicesDims(newaddr);
		backvec(dim) = backward(newaddr);
		dimvec(did) = dim;
		lev = indicesLevs(newaddr);
		newlevel(did) = lev;
		level(did) = lev;
		nsubgrids = nsubgrids * (uint32(lev)+1);
		down(did) = 1;
		newaddr = newaddr + 1;
		did = did + 1;
	end
	nsubgrids = nsubgrids - 1; 
	
	kend = kstart + npoints - 1;
	oldindex = currentindex;
	
	kl = uint32(1);
	while kl <= nsubgrids
		% Determine the subgrid to compute
		did = uint8(1);
		while did <= ndims
			if down(did)
				if level(did) > 0
					level(did) = level(did) - 1;
					oldindex = backvec(dimvec(did));
					break;
				else
					down(did) = 0;
				end
			else
				if level(did) < newlevel(did)
					level(did) = level(did) + 1;
					if directForward
						oldindex = forward(oldindex, dimvec(did));
					else						
						oldindex = forward(forwardAddr(oldindex), dimvec(did));
					end
					break;
				else
					down(did) = 1;
				end
			end
			did = did + 1;
		end
		
		oldndims = indicesNDims(oldindex);
		index = subGridAddr(oldindex);
		
		% Do special case
		if oldndims == 0
			k = kstart;
			while k <= kend
				ip(k) = ip(k) + z(index);
				k = k + 1;
			end
			kl = kl + 1;
			continue;	
		end
		
		oldaddr = indicesAddr(oldindex);
		noldpoints = subGridPoints(oldindex);
		
		% Find the matching old and new indices
		did = uint8(1);
		newdid = oldndims + 1;
		olddid = uint8(1);
		
		while did <= ndims
			if olddid <= oldndims
				olddim = indicesDims(oldaddr);
				if olddim == dimvec(did)
					% New and old dimension match
					olddimvec(olddid) = olddim;
					newdimvec(olddid) = olddim;
					oldlval = indicesLevs(oldaddr);
					oldlevel(olddid) = oldlval;
					newlevel2(olddid) = newlevel(did);
					backvec(olddim) = backward(oldaddr);
					oldaddr = oldaddr + 1;
					olddid = olddid + 1;
					did = did + 1;
					continue;
				end
			end
			% This part is reached if either the dimensions did not
			% match or there are no more old dimensions available	
			newdimvec(newdid) = dimvec(did);
			newlevel2(newdid) = newlevel(did);
			newdid = newdid + 1;
			did = did + 1;
		end
		
		[nord, order] = sort(oldlevel(1:oldndims),2,'descend');
		newdimvec(1:oldndims) = olddimvec(order);
		newlevel2(1:oldndims) = newlevel2(order);
		vals   = z(index:index+noldpoints-1);
	  ip(kstart:kend) = ip(kstart:kend) + ...
				spdctupstep(vals, nord, newlevel2(1:ndims), ...
										newdimvec(1:ndims));
		kl = kl + 1;
	end
	currentindex = currentindex + 1;
	kstart = kend + 1;
	nkl = nkl + 1;
end
