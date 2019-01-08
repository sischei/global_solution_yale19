function [ip,ipder,ipder2] = spcontdercc(d,z,y,levelseq,...
                               purgedata,maxlev)
% SPCONTDERCC   Compute interpolated values and continuous gradient 
%    vectors (Clenshaw-Curtis grid)
%    [IP,IPDER,IPRDER2] = SPCONTDERCC(D,Z,Y,LEVELSEQ,PURGEDATA,
%		 MAXLEV)  Computes interpolated values at grid values 
%    [Y1, ..., YN] and gradient vectors at augmented grid values 
%    to compute continuous derivatives via interpolation. 
%    (Internal function)
%
% See also SPDERIVCC

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : July 3, 2006

% Change log:
% V1.0   : July 3, 2006
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
index2 = zeros(d,1,'uint32');
index2d = zeros(d,2,'uint32'); 
repvec = ones(d+1,1,'uint32');
level  = ones(d,1,'uint8');
tempvec = zeros(d,1);
dervec  = zeros(d,2);
maxlev = double(maxlev);
stepsize = 1/2^maxlev;

% Conventional routine using levelseq as full array
% Get the number of levels
nlevels = uint32(size(levelseq,1));
	
if ~isempty(purgedata), purge = true; else purge = false; end
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
repvec(1) = 1;
for kl = 1:nlevels
	npoints = uint32(1);
	lval = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		level(k) = lval;
		if lval == 0
			repvec(k+1) = 1;
		elseif lval < 3
			repvec(k+1) = 2;
		else
			repvec(k+1) = 2^uint32(lval-1);
		end
		npoints = npoints * repvec(k+1);
		if k > 1
			repvec(k+1) = repvec(k+1) * repvec(k);
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
	ytd1 = 0;
	ytd2 = 0;
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
			
			% Compute augmented y values for derivatives
			if maxlev == 1
				ytd1 = 0.25;
				ytd2 = 0.75;
			elseif maxlev > 1
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
					index2(l) = 1;
					tempvec(l) = 1;
				else
					xp = floor(yt * 2);
					if xp == 0
						tempvec(l) = 2 * (0.5 - yt);
					else
						tempvec(l) = 2 * (yt - 0.5);
					end
					index2(l) = xp;
				end
				xp = floor(ytd1 * 2);
				if xp == 0
					dervec(l,1) = -2;
				else
					dervec(l,1) = 2;
				end
				index2d(l,1) = xp;
				xp = floor(ytd2 * 2);
				if xp == 0
					dervec(l,2) = -2;
				else
					dervec(l,2) = 2;
				end
				index2d(l,2) = xp;
			elseif lval == 0
				index2(l) = 0;
				tempvec(l) = 1;
				dervec(l,1) = 0;
				dervec(l,2) = 0;
				index2d(l,1) = 0;
				index2d(l,2) = 0;
			else
				scale = 2^double(lval);
				if yt == 1
					index2(l) = scale / 2 - 1;
					tempvec(l) = 0;
				else
					xp = floor(yt * scale / 2);
					index2(l) = xp;
					tempvec(l) = 1 - scale * abs(yt - (xp*2+1)/scale );
				end
				xp = floor(ytd1 * scale / 2);
				index2d(l,1) = xp;
				if ytd1 - ( (xp*2+1)/scale ) > 0
					dervec(l,1) = -scale;
				else
				  dervec(l,1) = scale;
				end
				xp = floor(ytd2 * scale / 2);
				index2d(l,2) = xp;
				if ytd2 - ( (xp*2+1)/scale ) > 0
					dervec(l,2) = -scale;
				else
				  dervec(l,2) = scale;
				end				
			end
			l = l + 1;
		end
		
		index3 = index + index2(1);
		for l = 2:d
			index3 = index3 + repvec(l)*index2(l);
		end
		temp = tempvec(1);
		for l = 2:d
			temp = temp * tempvec(l);
		end
		ip(k) = ip(k) + temp*z(index3);
		
		% Compute derivatives (left and right augmented one)
		for l = 1:d
			index3 = index;
			temp = 1;
			for l2 = 1:d
				if l == l2
					index3 = index3 + repvec(l2)*index2d(l2,1);
					temp = temp * dervec(l2,1);
				else
					index3 = index3 + repvec(l2)*index2(l2);
					temp = temp * tempvec(l2);
				end
			end
			ipder(k,l) = ipder(k,l) + temp*z(index3);
		end
		for l = 1:d
			index3 = index;
			temp = 1;
			for l2 = 1:d
				if l == l2
					index3 = index3 + repvec(l2)*index2d(l2,2);
					temp = temp * dervec(l2,2);
				else
					index3 = index3 + repvec(l2)*index2(l2);
					temp = temp * tempvec(l2);
				end
			end
			ipder2(k,l) = ipder2(k,l) + temp*z(index3);
		end
	end
	index = index + npoints;
end
