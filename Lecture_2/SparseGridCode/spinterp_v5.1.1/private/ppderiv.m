function [ipder] = ppderiv(ipder, ipder2, maxlev, y)
% PPDERIV  Derivative post processing. Combine two computed 
%    derivatives to one via linear interpolation.
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : July 1, 2006

% Change log:
% V1.1   : October 7, 2007
%          Fixed derivative computation such that when there
%          is a sign change, the derivative value is zero at
%          the grid point (important for optimization).
% V1.0   : July 1, 2006
%          Initial version

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

d = size(y,2);
ninterp = uint32(size(y,1));

if length(maxlev) < d
  % Expand maxLevel to a vector in case we have a
	% conventional, non-adaptive sparse grid.
	maxlev = ones(d,1,'uint8') * maxlev;
end

for k = 1:ninterp
	for l = 1:d
	  if maxlev(l) == 0, continue; end
		yt = y(k,l);
		stepsize = 1/2^double(maxlev(l));
		halfstep = 0.5 * stepsize;
		if maxlev(l) == 1
			ytd1 = 0.25;
		elseif maxlev(l) > 1
			if yt <= halfstep
				ytd1 = halfstep;
			elseif yt >= 1 - halfstep;
				ytd1 = 1 - 1.5*stepsize;
			else																							                     
				ytd1 = halfstep + ...
				       floor( (yt - halfstep) / stepsize) * stepsize;
			end
		end
		ipd1 = ipder(k,l);
		ipd2 = ipder2(k,l);
		if ipd1 * ipd2 >= 0 
			ipder(k,l) = ipd1 + (ipd2 - ipd1) / stepsize * (yt - ytd1);
		elseif yt <= ytd1 + halfstep
			ipder(k,l) = ipd1 - ipd1 / halfstep * (yt - ytd1);
		else
			ipder(k,l) = ipd2 / halfstep * (yt - ytd1 - halfstep);
		end
	end
end
