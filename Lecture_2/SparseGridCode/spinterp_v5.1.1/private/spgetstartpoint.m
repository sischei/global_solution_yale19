function [xopt, fval] = spgetstartpoint(z, xbox, options)
% SPGETSPARTPOINT   Computes the best start point for optimization
%    [XOPT,FVAL] = SPGETSTARTPOINT(Z, XBOX, OPTIONS)  Computes
%    best start point(s) for optimization from the available
%    sparse grid points, and eventually previous results and the
%    corner points for the box XBOX.
%    Optionally, spgetstartpoint also sets random start points or
%    user-specified start points (use OPTIONS created with 
%    SPOPTIMSET to configure this).
%
%    See also SPOPTIMSET
	 
% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.2
% Date   : September 1, 2007

% Change log:
% v1.2 : September 1, 2007
%        Added calculation of center function value if no grid
%        point selected.
% V1.1 : November 4, 2006
%        Added user-defined start point specification.
% V1.0 : June 9, 2005
%      : Initial version.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

if nargin < 3
	options = [];
end

minimize = spoptimget(options, 'Minimize', 'on');
maximize = spoptimget(options, 'Maximize', 'off');

isminimize = 0;
xoptmin = [];
ymin = [];
if strcmpi(minimize, 'on')
	isminimize = 1;
	ymin = inf;
end
	
ismaximize = 0;
xoptmax = [];
ymax = [];
if strcmpi(maximize, 'on')
	ismaximize = 1;
	ymax = -inf;
end

d = z.d;
% In case that no range has been provided to spvals -> set it to
% [0,1]^d. 
if isempty(z.range)
	z.range = [zeros(d,1) ones(d,1)];
end

startPoint = spoptimget(options, 'StartPoint', 'best');

if strcmpi(startPoint, 'best')
	prevResult = spoptimget(options, 'PrevResult', []);
	if ~isempty(prevResult)
		if isminimize
			ymin = prevResult(d+1,1);
			xoptmin = prevResult(1:d,1)';
			if ismaximize
				ymax = prevResult(d+1,2);
				xoptmax = prevResult(1:d,2)';
			end
		else
			ymax = prevResult(d+1,1);
			xoptmax = prevResult(1:d,1)';
		end
	end
	
	% Compute the extrema over the sparse grid points for the given
	% range
	[xopt, fval] = spfindopt(z, xbox, isminimize, ismaximize);
	
	% Compare to previous results; update if applicable.
	k = 1;
	if isminimize
		if fval(1) < ymin
			ymin = fval(1);
			xoptmin = xopt(:,1)';
		end
		k = 2;
	end
	if ismaximize
		if fval(k) > ymax
			ymax = fval(k);
			xoptmax = xopt(:,k)';
		end
	end
	
	% Compute the corner points if according option is set
	testCorners = spoptimget(options, 'TestCorners', 'off');
	if strcmpi(testCorners, 'on')
		xcell = num2cell(xbox,2);
		[xcell{1:d}] = ndgrid(xcell{:});
		y = spinterp(z, xcell{1:d});
		[fval, id] = max(y(:));
		
		% Compare to previous results; update if applicable.
		if isminimize
			if fval < ymin
				for k = 1:d
					xoptmin(k) = xcell{k}(id);
				end
			end
		end
		if ismaximize
			if fval > ymax
				for k = 1:d
					xoptmax(k) = xcell{k}(id);
				end
			end
		end
	end
	
	if isminimize & isempty(xoptmin)
		% if no point lies within the search box, take the center of
		% the box as start point
		xoptmin = ((xbox(:,1) + xbox(:,2))./2)';
		ymin = spsurfun(xoptmin, z);
	end
	
	if ismaximize & isempty(xoptmax)
		% if no point lies within the search box, take the center of
		% the box as start point
		xoptmax = ((xbox(:,1) + xbox(:,2))./2)';
		ymax = spsurfun(xoptmax, z);
	end
elseif strcmpi(startPoint, 'random')
	% Choose random start point in search box (may be shifted to a
  % full grid point later on).
	if isminimize
		xoptmin = (xbox(:,1) + rand(d,1).*(xbox(:,2)-xbox(:,1)))';
		ymin = spsurfun(xoptmin, z);
	end
	if ismaximize
		xoptmax = (xbox(:,1) + rand(d,1).*(xbox(:,2)-xbox(:,1)))';
		ymax = spsurfun(xoptmax, z);
	end
else
	% Use start point(s) provided by user
	k = 1;
	if size(startPoint,1) ~= d
		error('MATLAB:spinterp:badopt', ...
		  'Prodived start vector does not match sparse grid dimension');
	end
	if isminimize
		xoptmin = startPoint(:,k)';
		ymin = spsurfun(xoptmin, z);
	end
	if isminimize && ismaximize && size(startPoint,2) > 1
		k = k + 1;	
  end
	if ismaximize
		xoptmax = startPoint(:,k)';
		ymax = spsurfun(xoptmax, z);
	end
end
	
xopt = [xoptmin' xoptmax'];
fval = [ymin ymax];
