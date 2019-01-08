function x = getgpbaryw(allnx)
% Generate Gauss-Patterson barycentric weights for barypdstepgp
%    (internal function)

% Author : Andreas Klimke
% Version: 1.0
% Date   : November 18, 2007

% Change log:
% V1.0   : November 18, 2007
%					Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

d = uint8(length(allnx));
x = zeros(sum(allnx),1);

aid = uint32(0);
for k = 1:d
  x(1+aid:aid+allnx(k)) = gpbaryw(log2(double(allnx(k))+1)-1);
	aid = aid + allnx(k);
end
