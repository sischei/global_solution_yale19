function [ip, ipder] = spderivcc(d,z,y,levelseq,purgedata)
% SPDERIVCC   Compute interpolated values and gradient vectors 
%    (Clenshaw-Curtis grid)
%    [IP,IPDER] = SPDERIVCC(D,Z,Y,LEVELSEQ)  Computes interpolated
%    values and the gradients at [Y1, ..., YN] over the sparse grid 
%    for the given sequence of index sets LEVELSEQ. Y may be a 
%    double array for vectorized processing, with each row 
%    representing one point. 
%    The sparse grid data must be given as an array Z
%    containing the hierarchical surpluses (computed with
%    SPVALS). Note that the sparse grid is always normalized to the
%    unit cube [0,1]^D, i.e., if the weights have been computed for
%    a different domain, the values Y have to be rescaled
%    accordingly. (Internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : June 22, 2006

% Change log:
% V1.0   : June 22, 2006
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
index2 = zeros(d,1,'uint32'); 
repvec = ones(d,1,'uint32');
level  = ones(d,1,'uint8');
tempvec = zeros(d,1);
dervec  = zeros(d,1);

% Conventional routine using levelseq as full array
% Get the number of levels
nlevels = uint32(size(levelseq,1));
	
if ~isempty(purgedata), purge = true; else purge = false; end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
for kl = 1:nlevels
	npoints = uint32(1);
	lval = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		level(k) = lval;
		if lval == 0
			repvec(k) = 1;
		elseif lval < 3
			repvec(k) = 2;
		else
			repvec(k) = 2^uint32(lval-1);
		end
		npoints = npoints * repvec(k);
		if k > 1
			repvec(k) = repvec(k) * repvec(k-1);
		end
	end
	
	% Skip subgrids with all surpluses below droptol.
	if purge
		if purgedata(kl) == 0
			index = index + npoints;
			continue;
		end
	end
	
	yt = 0;
	for k = 1:ninterp
		temp = 1;
		dist = 0;
		l = uint16(1);
		l2 = uint16(1);
		while l <= d
			lval = level(l);
			yt = y(k,l);
			
			% security test if yt is within a valid range.
			if yt < 0  || yt > 1
				if yt < 0 yt = 0; else yt = 1; end
				warning('MATLAB:spinterp:outOfRange', ...
								'Interpolated point is out of valid range.');
			end
			
			% Compute the scaling factor and the array position of
			% the weight
			if lval == 1
				if yt == 1
					index2(l) = 1;
					tempvec(l) = 1;
					dervec(l) = 2;
				else
					xp = floor(yt * 2);
					if xp == 0
						tempvec(l) = 2 * (0.5 - yt);
						dervec(l) = -2;
					else
						tempvec(l) = 2 * (yt - 0.5);
						dervec(l) = 2;
					end
					index2(l) = xp;
				end
			elseif lval == 0
				index2(l) = 0;
				tempvec(l) = 1;
				dervec(l) = 0;
			else
				scale = 2^double(lval);
				if yt == 1
					index2(l) = scale / 2 - 1;
					tempvec(l) = 0;
					dervec(l) = -scale;
				else
					xp = floor(yt * scale / 2);
					index2(l) = xp;
					dist = yt - ( (xp*2+1)/scale );
					tempvec(l) = 1 - scale * abs(dist);
					
					% Note: It is important to have >= 0 and not something
					% like sign(dist) or dist > 0. It must be made sure that
					% the same left or right-sided derivative is taken for
					% all levels and bounds of basis functions to avoid
					% large jumps at the points of nondifferentiability.
					if dist >= 0 
						dervec(l) = -scale;
					else
						dervec(l) = scale;
					end
				end
			end
			l = l + 1;
		end
		
		index3 = index + index2(1);
		for l = 2:d
			index3 = index3 + repvec(l-1)*index2(l);
		end
		temp = tempvec(1);
		for l = 2:d
			temp = temp * tempvec(l);
		end
		ip(k) = ip(k) + temp*z(index3);
		
		% Compute derivatives
		for l = 1:d
			temp = 1;
			for l2 = 1:d
				if l == l2
					temp = temp * dervec(l2);
				else
					temp = temp * tempvec(l2);
				end
			end
			ipder(k,l) = ipder(k,l) + temp*z(index3);
		end
	end
	index = index + npoints;
end
