function [ip, ipder] = spderivcbsp(d,z,y,seq,purgedata)
% SPDERIVCBSP   Compute interpolated values and gradient vectors 
%    (Chebyshev grid, sparse indices)
%    [IP,IPDER] = SPDERIVCBSP(N,D,Z,Y)  Computes interpolated
%    values and the gradients at [Y1, ..., YN] over the sparse grid 
%    for the given sequence of index sets LEVELSEQ. Y may be a 
%    double array for vectorized processing, with each row 
%    representing one point. 
%    This routine does essentially the same as SPDERIVCB, but 
%    operates over sparse index sets for increased efficiency, 
%    especially in case of higher dimensions D. 
%    (Internal function)
%
% See also SPDERIVPCB.

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : June 30, 2006

% Change log:
% V1.0   : June, 30, 2006
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
	
if ~isempty(purgedata), purge = true; else purge = false; end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
index2 = zeros(d,1,'uint32'); 
level  = ones(d,1,'uint8');
dimvec = zeros(d,1,'uint16');
n = ones(d,1,'uint32');
np = ones(1,d,'uint32');
	
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
		n(did) = 2^uint32(lval)+1;
		if lval < 3
			np(did) = 2;
		else
			np(did) = 2^uint32(lval-1);
		end 
		npoints = subGridPoints(currentindex);
		addr = addr + 1;
		did = did + 1;
	end

	[t, order] = sort(level(1:ndims),1,'descend');
	nord = n(order);
	npord = np(order);
  order = uint16(dimvec(order));
	vals = z(index:index+npoints-1);
	ip = ip + barypdstepcb(vals, nord, order, getchebnodes(nord), y);
	
	% Compute derivatives
	if ndims > 1
		% make vals array ndims-dimensional
		valsnd = reshape(vals,npord);
		s.type = '()';
		s.subs(1:ndims) = {':'};
	end
	l = uint16(1);
	while l <= ndims
		if ndims > 1
			tmpvals = zeros(npord(l),ninterp);
			dimvec2 = [1:double(l)-1 double(l)+1:double(ndims)];
			tmpnord = nord(dimvec2);
			tmporder = order(dimvec2);
			for l2 = 1:npord(l)
				s.subs{l} = l2;
				vals = subsref(valsnd,s);
				tmpvals(l2,:) = barypdstepcb(vals(:), tmpnord, tmporder, ...
				                getchebnodes(tmpnord), y);
			end
			s.subs{l} = ':';
		else
			tmpvals = vals;
		end
		ipder(:,order(l)) = ipder(:,order(l)) + ...
		                    dctdiffcheb(nord(l), tmpvals, y(:,order(l)));
		l = l + 1;
	end
	
	index = index + npoints;
	currentindex = currentindex + 1;
end
