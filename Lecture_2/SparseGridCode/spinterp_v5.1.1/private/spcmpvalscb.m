function ip = spcmpvalscb(d,z,y,newlevelseq,levelseq)
% SPCMPVALSCB   Compute hierarchical surpluses, Chebyshev-grid
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.1
% Date   : January 24, 2006
	
% Change log:
% V1.0   : July 6, 2004
%          Initial release.
% V1.1   : January 24, 2006
%          Changed data types to operate on uint arrays. Corrected
%          Bug with nnewpoints array definition error, see
%          spcmpvalscc.m

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
n = ones(d,1,'uint32');
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
		if lval == 0
			n(k) = 1;
		else
			ndims = ndims + 1;
			n(k) = 2^uint32(lval)+1;
			if lval < 3
				npoints = npoints * 2;
			else
				npoints = npoints * 2^uint32(lval-1);
			end
		end
	end
  if npoints > 1
    [t, order] = sort(level,1,'descend');
    order = uint16(order(1:ndims));
    nord = n(order);
    x = getchebnodes(nord);
	  k = uint32(0);
	  for nkl = 1:nnewlevels
		  skiplevel = 0;
		  for l = 1:d
			  if level(l) > newlevelseq(nkl,l)
				  skiplevel = 1;
				  break;
			  end
		  end
		  kend = k + nnewpoints(nkl);
		
		  if ~skiplevel
			  vals = z(index:index+npoints-1);
			  ip(k+1:kend) = ip(k+1:kend) + ...
          barypdstepcb(vals, nord, order, x, y(k+1:kend,:));
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