function [ip,ipder,ipder2] = spcontderccsp(d,z,y,seq,...
                               purgedata,maxlevvec)
% SPCONTDERCCSP   Compute interpolated values and continuous 
%    gradient vectors (Clenshaw-Curtis grid)
%    [IP,IPDER,IPRDER2] = SPCONTDERCCSP(D,Z,Y,LEVELSEQ,PURGEDATA,
%    MAXLEVVEC)  Computes interpolated values at grid values 
%    [Y1, ..., YN] and gradient vectors at augmented grid values 
%    to compute continuous derivatives via interpolation. 
%    (Internal function)
%
% See also SPDERIVCCSP.

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : July 8, 2006

% Change log:
% V1.0   : July 8, 2006
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
ipder2 = zeros(ninterp,d);

indicesNDims  = seq.indicesNDims;
indicesDims   = seq.indicesDims;
indicesLevs   = seq.indicesLevs;
subGridPoints = seq.subGridPoints;
% Get the number of levels
nlevels = uint32(length(indicesNDims));
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);

if length(maxlevvec) < d
	maxlevvec = ones(d,1,'uint8') * maxlevvec;
end
	
if ~isempty(purgedata), purge = true; else purge = false; end
	
index2 = zeros(d,1,'uint32');
index2d = zeros(d,2,'uint32'); 
repvec = ones(d+1,1,'uint32');
level  = ones(d,1,'uint8');
dimvec = zeros(d,1,'uint16');
tempvec = zeros(d,1);
dervec  = zeros(d,2);
	
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
	repvec(1) = 1;
	while did <= ndims
		dimvec(did) = indicesDims(addr);
		lval = indicesLevs(addr);
		level(did) = lval;
		npoints = subGridPoints(currentindex);
		if lval < 3
			repvec(did+1) = 2;
		else
			repvec(did+1) = 2^(uint32(lval)-1);
		end
		if did > 1
			repvec(did+1) = repvec(did+1) * repvec(did);
		end
		addr = addr + 1;
		did = did + 1;
	end
	
	yt = 0;
	ytd1 = 0;
	ytd2 = 0;
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
				
			% Compute augmented y values for derivatives
			maxlev = double(maxlevvec(l));
			if maxlev == 1
				ytd1 = 0.25;
				ytd2 = 0.75;
			elseif maxlev > 1
				stepsize = 1/2.^maxlev;
				halfstep = 0.5 * stepsize;
				if yt <= halfstep
					ytd1 = halfstep;
					ytd2 = 1.5 * stepsize;
				elseif yt >= 1 - halfstep;
					ytd1 = 1 - 1.5*stepsize;
					ytd2 = 1 - halfstep;
				else																							                     
					ytd1 = halfstep + ...
					        floor( (yt - halfstep) / stepsize) * stepsize;
					ytd2 = ytd1 + stepsize;
				end
			end
			
			% Compute the scaling factor and the array position of
			% the weight
			if lval == 1
				if yt == 1
					index2(did) = 1;
					tempvec(did) = 1;
				else
					xp = floor(yt * 2);
					if xp == 0
						tempvec(did) = 2 * (0.5 - yt);
					else
						tempvec(did) = 2 * (yt - 0.5);
					end
					index2(did) = xp;
				end
				xp = floor(ytd1 * 2);
				if xp == 0
					dervec(did,1) = -2;
				else
					dervec(did,1) = 2;
				end
				index2d(did,1) = xp;
				xp = floor(ytd2 * 2);
				if xp == 0
					dervec(did,2) = -2;
				else
					dervec(did,2) = 2;
				end
				index2d(did,2) = xp;
			else
				scale = 2^double(lval);
				if yt == 1
					index2(did) = scale / 2 - 1;
					tempvec(did) = 0;
				else
					xp = floor(yt * scale / 2);
					index2(did) = xp;
					tempvec(did) = 1 - scale * abs(yt - (xp*2+1)/scale );
				end
				xp = floor(ytd1 * scale / 2);
				index2d(did,1) = xp;
				if ytd1 - ( (xp*2+1)/scale ) > 0
					dervec(did,1) = -scale;
				else
				  dervec(did,1) = scale;
				end
				xp = floor(ytd2 * scale / 2);
				index2d(did,2) = xp;
				if ytd2 - ( (xp*2+1)/scale ) > 0
					dervec(did,2) = -scale;
				else
				  dervec(did,2) = scale;
				end				
			end
			did = did + 1;
		end
			
		index3 = index + index2(1);
		did = uint8(2);
		while did <= ndims
			index3 = index3 + repvec(did)*index2(did);
			did = did + 1;
		end
		temp = tempvec(1);
		did = uint8(2);
		while did <= ndims
			temp = temp * tempvec(did);
			did = did + 1;
		end
		ip(k) = ip(k) + temp*z(index3);

		% Compute derivatives (left and right augmented one)
		did = uint8(1);
		while did <= ndims
			index3 = index;
			did2 = uint8(1);
			temp = 1;
			while did2 <= ndims
				if did == did2
					index3 = index3 + repvec(did2)*index2d(did2,1);
					temp = temp * dervec(did2,1);
				else
					index3 = index3 + repvec(did2)*index2(did2);
					temp = temp * tempvec(did2);
				end
				did2 = did2 + 1;
			end
			l = dimvec(did);
			ipder(k,l) = ipder(k,l) + temp*z(index3);
			did = did + 1;
		end
		did = uint8(1);
		while did <= ndims
			index3 = index;
			did2 = uint8(1);
			temp = 1;
			while did2 <= ndims
				if did == did2
					index3 = index3 + repvec(did2)*index2d(did2,2);
					temp = temp * dervec(did2,2);
				else
					index3 = index3 + repvec(did2)*index2(did2);
					temp = temp * tempvec(did2);
				end
				did2 = did2 + 1;
			end
			l = dimvec(did);
			ipder2(k,l) = ipder2(k,l) + temp*z(index3);
			did = did + 1;
		end

		k = k + 1;
	end
	index = index + npoints;
	currentindex = currentindex + 1;
end
