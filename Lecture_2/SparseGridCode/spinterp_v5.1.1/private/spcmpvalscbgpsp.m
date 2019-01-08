function ip = spcmpvalscbgpsp(d,z,y,seq,fromindex,toindex,isgp)
% SPCMPVALSCBGPSP   Compute surpluses, Chebyshev or Gauss-
%    Patterson grid, sparse indices (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.2
% Date   : December 10, 2007
	
% Change log:
% V1.0   : March 4, 2005 
%          Initial release.
% V1.1   : March 28, 2005
%          Added indirect addressing of forward neighbor arrays.
% v1.2   : December 10, 2007
%          Created from spcmpvalscbsp and added capability of
%          handling Gauss-Patterson grid type.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

% Check if isgp is provided
if nargin < 7, isgp = false; end

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
index2 = zeros(d,1,'uint32'); 
level = zeros(d,1,'uint8');
newlevel = zeros(d,1,'uint8');
oldlevel = zeros(d,1,'uint8');
dimvec = ones(d,1,'uint16');
olddimvec = ones(d,1,'uint16');
down = ones(d,1,'uint8');
n    = ones(d,1,'uint32');

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
		did = uint16(1);
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
		did = uint16(1);
		while did <= oldndims
			olddim = indicesDims(oldaddr);
			olddimvec(did) = olddim;
			oldlval = indicesLevs(oldaddr);
			oldlevel(did) = oldlval;
			backvec(olddim) = backward(oldaddr);
			if isgp
			  n(did) = 2^uint32(oldlval+1) - 1;
			else
				n(did) = 2^uint32(oldlval) + 1;
			end
			oldaddr = oldaddr + 1;
			did = did + 1;
		end
	  vals  = z(index:index+noldpoints-1);
		if isgp
		  nord  = n(1:oldndims);
			order = uint16(olddimvec(1:oldndims));
			ip(kstart:kend) = ...
  			ip(kstart:kend) + barypdstepgp(vals, nord, order, ...
        getgpnodes(nord), y(kstart:kend,:), getgpbaryw(nord));
		else
			[t, order] = sort(oldlevel(1:oldndims),1,'descend');
			nord   = n(order);
			order  = uint16(olddimvec(order));
			ip(kstart:kend) = ...
				ip(kstart:kend) + barypdstepcb(vals, nord, order, ...
        getchebnodes(nord), y(kstart:kend,:));
		end
		kl = kl + 1;
	end
	currentindex = currentindex + 1;
	kstart = kend + 1;
	nkl = nkl + 1;
end
