function [isdispiter, iterstr] = initoptidisp(dispopt)
% INITOPTIDISP   Initializes iteration information display
%   (Internal function)

% Author : Andreas Klimke
% Date   : October 28, 2006
% Version: 1.0

% Change log:
% V1.0   : October 28, 2006
%          Initial version

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

if strcmpi(dispopt, 'iter') 
	disp(' Iteration   Func-count Grad-count     f(x)            Procedure');
	iterstr =' %5.0f        %5.0f     %5.0f     %12.6g         %s';
	isdispiter = true;
else
	iterstr = [];
	isdispiter = false;
end
