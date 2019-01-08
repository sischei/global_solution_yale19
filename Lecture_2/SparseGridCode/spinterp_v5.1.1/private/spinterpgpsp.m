function ip = spinterpgpsp(d,z,y,seq,purgedata)
% SPINTERPGPSP  Polynomial interpolation (Gauss-P., sparse indices)
%    IP = SPINTERPGP(D,Z,Y,SEQ)  Computes the interpolated
%    values at [Y1, ..., YN] over the sparse grid for the given
%    sequence of index sets SEQ. at level N. This routine does 
%    essentially the same as SPINTERPGP, but operates over sparse 
%    index sets for increased efficiency, especially in case of 
%    higher dimensions D. (Internal function)
%
% See also SPINTERPGP.

% Author : Andreas Klimke
% Version: 1.0
% Date   : December 22, 2007

% Change log:
% V1.0   : December 22, 2007
%          Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

ninterp = uint32(size(y,1));
ip = zeros(ninterp,1);

indicesNDims  = seq.indicesNDims;
indicesDims   = seq.indicesDims;
indicesLevs   = seq.indicesLevs;
subGridPoints = seq.subGridPoints;
% Get the number of levels
nlevels = uint32(length(indicesNDims));
	
if ~isempty(purgedata), purge = true; else purge = false; end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
index2 = zeros(d,1,'uint32'); 
level  = ones(d,1,'uint8');
dimvec = zeros(d,1,'uint16');
n = ones(d,1,'uint32');
	
% currentindex is the index of the currently processed multiindex
% entry. 
currentindex = uint32(1);
addr = uint32(1);

while currentindex <= nlevels
	npoints = 1;
	lval = 0;
	ndims = indicesNDims(currentindex);
	% Do special case
	if ndims == 0
		k = uint32(1);
		while k <= ninterp
			ip(k) = ip(k) + z(index);
			k = k + 1;
		end
		index = index + 1;
		currentindex = currentindex + 1;
		continue;
	end
	
	% Skip subgrids with all surpluses below droptol.
	if purge
		if purgedata(currentindex) == 0
			index = index + subGridPoints(currentindex);
			addr = addr + uint32(ndims);
			currentindex = currentindex + 1;
			continue;
		end
	end
	
	npoints = subGridPoints(currentindex);
	did = uint8(1);
	while did <= ndims
		dimvec(did) = indicesDims(addr);
		lval = indicesLevs(addr);
		level(did) = lval;
		n(did) = 2^uint32(lval+1)-1;
		addr = addr + 1;
		did = did + 1;
	end

	nord = n(1:ndims);
  order = uint16(dimvec(1:ndims));
	vals = z(index:index+npoints-1);
	ip = ip + barypdstepgp(vals, nord, order, getgpnodes(nord), y, ...
	  getgpbaryw(nord));
	
	index = index + npoints;
	currentindex = currentindex + 1;
end
