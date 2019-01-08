function ip = spcmpvalsccsp(d,z,y,seq,fromindex,toindex)
% SPCMPVALSCCSP   Compute surpluses, Clenshaw-C. grid, sparse indices
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.1
% Date   : March 28, 2005
	
% Change log:
% V1.0   : March 4, 2005 
%          Initial release.
% V1.1   : March 28, 2005
%          Added indirect addressing of forward neighbor arrays.

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
index2 = zeros(d,1,'uint32'); 
repvec = ones(d,1,'uint32');
level = zeros(d,1,'uint8');
backvec = ones(d,1,'uint32');
newlevel = zeros(d,1,'uint8');
oldlevel = zeros(d,1,'uint8');
dimvec = ones(d,1,'uint16');
olddimvec = ones(d,1,'uint16');
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
			if oldlval < 3
				repvec(did) = 2;
			else
				repvec(did) = 2^(uint32(oldlval)-1);
			end
			if did > 1
				repvec(did) = repvec(did) * repvec(did-1);
			end
			oldaddr = oldaddr + 1;
			did = did + 1;
		end

		k = kstart;
		while k <= kend
			temp = 1;
			did  = uint16(1); 
			while did <= oldndims
				l = olddimvec(did);
				lval = oldlevel(did);
				yt = y(k,l);
				
				if lval == 1
					if yt == 1
						index2(did) = 1;
					else
						xp = floor(yt * 2);
						if xp == 0
							temp = temp * 2 * (0.5 - yt);
						else
							temp = temp * 2 * (yt - 0.5);
						end
						index2(did) = xp;
					end
				else
					scale = 2^double(lval);
					xp = floor(yt * scale / 2);
					temp = temp * ...
								 (1 - scale * abs( yt - ( (xp*2+1)/scale )));
					index2(did) = xp;
				end
				did = did + 1;
			end
					
			%	if temp > 0, etc, has been removed, since in case of
			%	computing the hierarchical surpluses, when omitting the
			% non-contributing subdomands, this case cannot occur anyway.
			index3 = index + index2(1);
			did = uint16(2);
			while did <= oldndims
				index3 = index3 + repvec(did-1)*index2(did);
				did = did + 1;
			end
			ip(k) = ip(k) + temp*z(index3);
			k = k + 1;
		end
		kl = kl + 1;
	end
	currentindex = currentindex + 1;
	kstart = kend + 1;
	nkl = nkl + 1;
end
	