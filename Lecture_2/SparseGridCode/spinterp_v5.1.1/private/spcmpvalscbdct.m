function ip = spcmpvalscbdct(d,z,y,newlevelseq,levelseq)
% SPCMPVALSCBFFT   Compute surpluses using DCT. Otherwise,
%    does the same as SPCMPVALSCB
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : February 10, 2006
	
% Change log:
% V1.0   : February 10, 2006
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

% Get the number of new levels
nnewlevels = uint32(size(newlevelseq, 1));
nnewpoints = zeros(nnewlevels,1,'uint32');

% Get the number of old levels
nlevels = uint32(size(levelseq,1));
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = uint32(1);
	
index2 = zeros(d,1,'uint32'); 
n = ones(1,d,'uint32');
nnew = ones(1,d,'uint32');
dims = zeros(1,d,'uint16');
newdims = zeros(1,d,'uint16');
level = ones(d,1,'uint8');

% Compute the number of points per subdomain of new levels
for kl = 1:nnewlevels
	npoints = uint32(1);
	lval = uint8(0);
	for k = 1:d
		lval = newlevelseq(kl,k);
		if lval == 0 % do nothing
		elseif lval < 3
			npoints = npoints * 2;
		else
			npoints = npoints * 2^uint32(lval-1);
		end
	end
	nnewpoints(kl) = npoints;
end	

for kl = 1:nlevels
	npoints = uint32(1);
	lval = uint8(0);
	ndims = uint8(0);
	for k = 1:d
		lval = levelseq(kl,k);
		level(k) = lval;
		if lval > 0
			ndims = ndims + 1;
			dims(ndims) = k;
			n(ndims) = lval;
			if lval < 3
				npoints = npoints * 2;
			else
				npoints = npoints * 2^uint32(lval-1);
			end
		end
	end
	nolddims = ndims;
  if npoints > 1
    [t, order] = sort(n(1:ndims),2,'descend');
    n(1:ndims) = n(order);
    dims(1:ndims) = dims(order);
	  k = uint32(0);
	  for nkl = 1:nnewlevels
		  skiplevel = 0;
		  ndims = nolddims;
		  for l = 1:d
				lval = level(l);
				newlval = newlevelseq(nkl,l);
			  if lval > newlval
				  skiplevel = 1;
				  break;
			  elseif lval < newlval
					nnew(l) = newlval;
					if lval == 0
						ndims = ndims + 1;
						dims(ndims) = l;
					end
				else
					nnew(l) = newlval;
				end
		  end
		  kend = k + nnewpoints(nkl);
		  if ~skiplevel
			  vals = z(index:index+npoints-1);
			  ip(k+1:kend) = ip(k+1:kend) + ...
						spdctupstep(vals, n(1:nolddims), nnew(dims(1:ndims)), ...
												dims(1:ndims));
			end
		  k = kend;
	  end
  else
    % no interpolation necessary; we have just a single point with a 
    % constant function
    ip = ip + z(index);
  end
	index = index + npoints;
end
