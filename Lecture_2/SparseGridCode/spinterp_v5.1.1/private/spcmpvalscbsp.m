function ip = spcmpvalscbsp(d,z,y,seq,fromindex,toindex)
% SPCMPVALSCBSP   Compute surpluses, Chebyshev grid, sparse indices
%    (internal function)

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : December 10, 2007
	
% Change log:
% V1.0   : December 10, 2007 
%          Wrapper to spcmpvalscbgpsp for Chebyshev call.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

ip = spcmpvalscbgpsp(d,z,y,seq,fromindex,toindex);
