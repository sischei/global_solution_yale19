function seq = spgetseqsp(n, d, options)
% SPGETSEQSP   Get the sequence of index sets in sparse format
%    (internal function)

% Author : Andreas Klimke
% Date   : November 18, 2007
% Version: 1.1

% V1.1   : November 18, 2007
%          Added new grid type : Gauss-Patterson

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

n = uint8(n);
d = uint16(d);

if nargin < 3, options = []; end

prevseq = spget(options, 'PrevResults');
gridtype = spget(options, 'GridType', 'Clenshaw-Curtis');

if strcmpi(gridtype, 'maximum')
	error('MATLAB:spinterp:badopt',['Sparse indices are not supported' ...
										' with maximum-norm based grid.']);
end

if isempty(prevseq)
	% Total of levels to compute
	nlevels    = uint32(nchoosek(double(n)+double(d),double(d)));
	if n > 0
		noldlevels = uint32(nchoosek(double(n)+double(d)-1,double(d)));
	else
		noldlevels = uint32(0);
	end
	min_nd = min(double(n),double(d));
	
	% Indices data
	indicesNDims = zeros(nlevels,1,'uint8');
	indicesDims  = zeros(nlevels*min_nd,1,'uint16');
	indicesLevs  = zeros(nlevels*min_nd,1,'uint8');
	indicesAddr  = zeros(nlevels,1,'uint32');
	
	% Neighbor data
	backward     = zeros(nlevels*min_nd,1,'uint32');
	forward      = zeros(noldlevels,d,'uint32');
	
	% Indicate active/passive indices (only necessary for
  % construction of the grid)
	active    = false(nlevels,1);
	
	% Store start adresses of a sub-grid in the support node array
	ptsaddr   = zeros(nlevels,1, 'uint32');	
	ptslen    = zeros(nlevels,1, 'uint32');
	
	% Initialize with first index; There are no values to be inserted
  % in the indices array, since the first index is (0,...,0). (We
  % start from 0 -- this also results in the sparse structure of
  % the array).
	currentindex = uint32(1);
	ni = uint32(1);

	ptsaddr(1) = 1;
	indicesAddr(1) = 1;
	indicesDims(1) = 0;

	% Compute number of nodes of first sub-grid, and set grid code
	% for computation of other sub-grid sizes in loop.
	switch lower(gridtype)
	 case 'clenshaw-curtis'
		ptslen(1) = uint32(1);
	 case 'chebyshev'
		ptslen(1) = uint32(1);
	 case 'maximum'
		ptslen(1) = uint32(3^d);
	 case 'noboundary'
		ptslen(1) = uint32(1);
	 case 'gauss-patterson'
		ptslen(1) = uint32(1);
	end
	nold = uint32(0);
else
	nold = prevseq.maxN;
	
	% Total of new levels to compute
	nlevels = uint32(nchoosek(double(d)+double(n),double(d)) ...
									 - size(prevseq.indicesNDims,1));
	noldlevels = uint32(nchoosek(double(d)+double(n)-1, double(d)) ...
															 - size(prevseq.forwardNeighbors, ...
																			1));
	
	
	min_nd = min(double(n),double(d));

	% Indices data
	indicesNDims = [prevseq.indicesNDims; zeros(nlevels,1,'uint8')];
	indicesDims  = [prevseq.indicesDims; ...
									zeros(nlevels*min_nd,1,'uint16')];
	indicesLevs  = [prevseq.indicesLevs; ...
									zeros(nlevels*min_nd,1,'uint8')];
	indicesAddr  = [prevseq.indicesAddr; ...
									zeros(nlevels,1,'uint32')];
	
	% Neighbor data
	backward     = [prevseq.backwardNeighbors; ...
									zeros(nlevels*min_nd,1,'uint32')];
	forward      = [prevseq.forwardNeighbors; ...
									zeros(noldlevels,d,'uint32')];
	
	% Indicate active/passive indices (only necessary for
  % construction of the grid)
	active    = [prevseq.activeIndices; false(nlevels,1)];
	
	% Store start adresses of a sub-grid in the support node array
	ptsaddr   = [prevseq.subGridAddr; zeros(nlevels,1, 'uint32')];	
	ptslen    = [prevseq.subGridPoints; zeros(nlevels,1, 'uint32')];

	% Initialize with first index; There are no values to be inserted
  % in the indices array, since the first index is (0,...,0). (We
  % start from 0 -- this also results in the sparse structure of
  % the array).
	currentindex = prevseq.currentindex;
	ni = uint32(size(prevseq.indicesNDims,1));
end


% Compute number of nodes of first sub-grid, and set grid code
% for computation of other sub-grid sizes in loop.
switch lower(gridtype)
 case 'clenshaw-curtis'
	gridcode = 0;
 case 'chebyshev'
	gridcode = 0;
 case 'maximum'
	gridcode = 1;
 case 'noboundary'
	gridcode = 2;
 case 'gauss-patterson'
	gridcode = 2;
end

% Allocate some temporary arrays
bn    = zeros(d,1,'uint32');
bdim  = zeros(d,1,'uint16');
bid   = zeros(d,1,'uint32');
level = zeros(d,1,'uint8');

l = nold;
addr = uint32(1);

while l < n
	nnewlevels = uint32(nchoosek(double(d)+double(l)-1,double(d)-1));
	k = uint32(1);
	while k <= nnewlevels
		active(currentindex) = false;
		
		% Get all the backward neighbors
		nbid = indicesNDims(currentindex);
		addr = indicesAddr(currentindex);
		
		j = uint8(1);
		while j <= nbid
			bdim(j)  = indicesDims(addr);
			bid(j)   = backward(addr);
			level(bdim(j)) = indicesLevs(addr);
			addr = addr + 1;
			j = j + 1;
		end
		
		i = uint16(1);
		while i <= d
			
			isadmissible = true;
			j = uint8(1);
			while j <= nbid
				if i == bdim(j) 
					j = j + 1;
					continue; 
				end
				backofnew = forward(bid(j),i);
				if active(backofnew)
					isadmissible = false;
					% TODO: It seems that it is OK to set i = d here. Before:
          % i = i + 1 which was definitely wrong; but worked anyway.
					i = d;
					break;
				else
					bn(bdim(j)) = backofnew;
				end
				j = j + 1;
			end
			
			if isadmissible
				bn(i) = currentindex;

        addr = indicesAddr(ni)+uint32(indicesNDims(ni));
        ni = ni + 1;
				indicesAddr(ni) = addr;
				
				if nbid == 0
					nnewdims = uint8(1);
					indicesDims(addr) = i;
					indicesLevs(addr) = uint8(1);
					backward(addr) = bn(i);
					forward(bn(i),i) = ni;
					npoints = uint32(2);
					addr = addr + 1;
				else
					j = uint8(1);
					insert = true;
					nnewdims = uint8(0);
					npoints = uint32(1);
					while j <= nbid
						nnewdims = nnewdims + 1;
						id = bdim(j);
						if i > id
							lev = level(id);
							j  = j + 1;
						elseif i == id
							insert = false;
							lev = level(id) + 1;
							j = j + 1;
						elseif insert
							insert = false;
							lev = uint8(1);
							id = i;
						else
							lev = level(id);
							j = j + 1;
						end
						forward(bn(id),id) = ni;
						indicesDims(addr) = id;
						indicesLevs(addr) = lev;
						backward(addr) = bn(id);
						
						% Compute the number of gridpoints of according
            % sub-grid
						if gridcode == 0
							if lev < 3
								npoints = npoints * 2;
							else
								npoints = npoints * 2^uint32(lev-1);
							end
						elseif gridcode == 1
							if lev == 0
								npoints = npoints * 3;
							else
								npoints = npoints * 2^uint32(lev);
							end
						else
							npoints = npoints * 2^uint32(lev);
						end
						addr = addr + 1;
					end
				end
				indicesNDims(ni) = nnewdims;
				ptslen(ni) = npoints;
				ptsaddr(ni) = ptsaddr(ni-1) + ptslen(ni-1);		
				active(ni) = true;
			end
			i = i + 1;
		end
		currentindex = currentindex + 1;
		k = k + 1;
	end
	l = l + 1;
end

seq.indicesNDims = indicesNDims;
seq.indicesDims  = indicesDims(1:addr-1);
seq.indicesLevs  = indicesLevs(1:addr-1);
seq.indicesAddr  = indicesAddr;
seq.activeIndices = active; 
seq.forwardNeighbors = forward; 
seq.backwardNeighbors = backward(1:addr-1); 
seq.subGridPoints = ptslen; 
seq.subGridAddr = ptsaddr; 
seq.maxN = n; 
seq.currentindex = currentindex;
