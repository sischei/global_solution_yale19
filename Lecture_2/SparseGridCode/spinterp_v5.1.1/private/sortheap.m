function A = sortheap(A, na, G)
% SORTHEAP  Reorder heap after a new element has been added at end
%    (Internal function)
		
% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : April 20, 2004

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

k = na;
while k > 1
	prev = floor(k/2);
	if G(A(prev)) < G(A(k))
		temp = A(prev);
		A(prev) = A(k);
		A(k) = temp;
	else
		break;
	end
	k = floor(k/2);
end
