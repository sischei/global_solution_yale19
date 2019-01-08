function [ip, ipder] = spinterpccsp(d,z,y,seq,purgedata)
% SPDERIVPCCSP  Compute interpolated values and gradients
%    (Clenshaw-Curtis, sparse ids)
%    [IP,IPDER] = SPINTERPCC(D,Z,Y,SEQ)  Computes interpolated
%    values and the gradients at [Y1, ..., YN] over the sparse
%    grid for the given sequence of index sets SEQ. This routine
%    does essentially the same as SPDERIVCC, but operates over
%    sparse index sets for increased efficiency, especially in
%    case of higher dimensions D. (Internal function)
%
% See also SPDERIVCC.

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : June 23, 2006

% Change log:
% V1.0   : June 23, 2006
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
ipder = zeros(ninterp,d);

indicesNDims  = seq.indicesNDims;
indicesDims   = seq.indicesDims;
indicesLevs   = seq.indicesLevs;
subGridPoints = seq.subGridPoints;
% Get the number of levels
nlevels = uint32(length(indicesNDims));
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
if ~isempty(purgedata), purge = true; else purge = false; end
	
index2 = zeros(d,1,'uint32'); 
repvec = ones(d,1,'uint32');
level  = ones(d,1,'uint8');
dimvec = zeros(d,1,'uint16');
tempvec = zeros(d,1);
dervec  = zeros(d,1);
	
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
	
	did = uint8(1);
	while did <= ndims
		dimvec(did) = indicesDims(addr);
		lval = indicesLevs(addr);
		level(did) = lval;
		npoints = subGridPoints(currentindex);
		if lval < 3
			repvec(did) = 2;
		else
			repvec(did) = 2^(uint32(lval)-1);
		end
		if did > 1
			repvec(did) = repvec(did) * repvec(did-1);
		end
		addr = addr + 1;
		did = did + 1;
	end
	yt = 0;
		
	k = uint32(1);
	while k <= ninterp
		temp = 1;
		dist = 0;
		did = uint8(1);
		did2 = uint8(1);
		while did <= ndims
			l = dimvec(did);
			lval = level(did);
			yt = y(k,l);
			
			% security test if yt is within a valid range.
			if yt < 0  || yt > 1
				if yt < 0, yt = 0; else yt = 1; end
				warning('MATLAB:spinterp:outOfRange', ...
								'Interpolated point is out of valid range.');
			end
				
			% Compute the scaling factor and the array position of
			% the weight
			if lval == 1
				if yt == 1
					index2(did) = 1;
					tempvec(did) = 1;
					dervec(did) = 2;
				else
					xp = floor(yt * 2);
					if xp == 0
						tempvec(did) = 2 * (0.5 - yt);
						dervec(did) = -2;
					else
						tempvec(did) = 2 * (yt - 0.5);
						dervec(did) = 2;
					end
					index2(did) = xp;
				end
			else
				scale = 2^double(lval);
				if yt == 1
					index2(did) = scale / 2 - 1;
					tempvec(did) = 0;
					dervec(did) = -scale;
				else
					xp = floor(yt * scale / 2);
					index2(did) = xp;
					dist = yt - ( (xp*2+1)/scale );
					tempvec(did) = 1 - scale * abs(dist);
					if dist >= 0
						dervec(did) = -scale;
					else
						dervec(did) = scale;
					end
				end
			end
			did = did + 1;
		end
			
		index3 = index + index2(1);
		did = uint8(2);
		while did <= ndims
			index3 = index3 + repvec(did-1)*index2(did);
			did = did + 1;
		end
		temp = tempvec(1);
		did = uint8(2);
		while did <= ndims
			temp = temp * tempvec(did);
			did = did + 1;
		end
		ip(k) = ip(k) + temp*z(index3);
		
		% Compute derivatives
		did = uint8(1);
		while did <= ndims
			did2 = uint8(1);
			temp = 1;
			while did2 <= ndims
				if did == did2
					temp = temp * dervec(did2);
				else
					temp = temp * tempvec(did2);
				end
				did2 = did2 + 1;
			end
			l = dimvec(did);
			ipder(k,l) = ipder(k,l) + temp*z(index3);
			did = did + 1;
		end

		k = k + 1;
	end
	index = index + npoints;
	currentindex = currentindex + 1;
end
