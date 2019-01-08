function [index, A] = popheap(A, na, G)
% POPHEAP  Pop first element out of heap
%    (internal function)
		
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

index = A(1);

% reorder the remaining heap
k = 1;
while k*2 < na
	next = k * 2;
	if G(A(next)) > G(A(next+1))
		A(k) = A(next);
		k = next;
	else
		A(k) = A(next + 1);
		k = next + 1;
	end
end
% move the last element of A to the free spot
A(k) = A(na);

% move that element up the heap to the correct position.
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
