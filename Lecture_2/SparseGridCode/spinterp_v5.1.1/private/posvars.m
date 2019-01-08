function val = posvars(d, varpos, varargin)
% POSVARS   Shift interpolation parameters to right position
%    Shifts the interpolation parameters to the right position, and 
%    fills field up with any parameters in varargin. (Internal 
%    function).	

% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.0
% Date   : April 10, 2004

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

if ~isempty(varpos)
	val = cell(1,d+length(varargin));
	for l = 1:d
		val{varpos(l)} = 1;
	end
	m = 1;
	if length(varargin) > 0
		for l = 1:d+length(varargin)
			if isempty(val{l})
				val{l} = varargin{m};
				m = m + 1;
			end
		end
	end
else
	val = varargin;
end
