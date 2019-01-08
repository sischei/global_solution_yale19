function x = getchebnodes(allnx)
% Generate Chebyshev nodes for barypdstepcb
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : January 25, 2006

% Change log:
% V1.0   : Initial release.

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
  x(1+aid:aid+allnx(k)) = 0.5 - cos(linspace(0,1,allnx(k)) * pi) / 2;
	aid = aid + allnx(k);
end
