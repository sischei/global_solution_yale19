function [xopt, fval] = spfindopt(z, xbox, isminimize, ismaximize)
% SPFINTOPT   Compute the extrema over the sparse grid points for
% the box XBOX
	 
% Author : Andreas Klimke, Universitaet Stuttgart
% Version: 1.5
% Date   : May 14, 2008

% Change log:
% V1.0 : February 18, 2004
%      : Initial version.
% V1.1 : September 6, 2004
%        Added capability for handling dimension-adaptive data.
% V1.2 : June 9, 2005
%        Added z.indices.currentindex = 1 statement to handle
%        sparse indices arrays correctly.
% V1.3 : August 1, 2005
%        Corrected bug; Stored function values were not used
%        if grid was not stored due to wrong placement of if
%        statement.
% V1.4 : September 1, 2007
%        Fixed bug in cropping (used and instead of or).
% V1.5 : May 14, 2008
%        Fixed bug in cropping (added tolerance, fixed indices
%        range).

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

if nargin < 3;
  isminimize = 1;
  ismaximize = 0;
end

if isfield(z, 'indices')
	sparseIndices = 1;
else
	sparseIndices = 0;
end

if isfield(z, 'selectOutput')
	output = z.selectOutput;
else
	output = 1;
end

xoptmin = []; xoptmax = [];
ymin = []; ymax = [];

if isminimize, ymin = inf; end
if ismaximize, ymax = -inf; end

d = z.d;

n = size(z.vals, 2);

for k = 1:n
	% Get the sparse grid points
	y = [];
	if isfield(z, 'fvals')
		y = z.fvals{output, k};
	end	
	if isfield(z, 'grid')
		sx = z.grid{k};
	else
		if ~sparseIndices
			sx = spgrid(k-1,d,spset('GridType', z.gridType));
		else
			% Make sure that the entire grid is returned!
			% TODO: This should be broken down in the future for very
      % high-dimensional interpolants and many support nodes, since
      % it may require to much memory do store the entire grid in
      % one single array.
			z.indices.currentindex = 1;
			sx = spgrid(z.indices(k,:),[],spset('GridType', z.gridType));
		end
		for l = 1:d
			% Rescale sparse grid to actual range
			sx(:,l) = z.range(l,1) + (z.range(l,2)-z.range(l,1)).*sx(:,l);
		end
	end
	% Crop all points outside of search box
	skip = 0;
	np = size(sx,1);
	id = zeros(np,1);
	nid = 0;
	for k = 1:np
	  crop = 0;
		for l = 1:d
			if sx(k,l) < xbox(l,1) - 10*eps*(z.range(l,2)-z.range(l,1)) ...
			   || sx(k,l) > xbox(l,2) + 10*eps*(z.range(l,2)-z.range(l,1))
			  crop = 1;
				break;
			end
		end
		if crop == 0
			nid = nid + 1;
			id(nid) = k;
		end
	end
  sx = sx(id(1:nid),:);
  
	if ~isempty(sx)
		if ~isempty(y)
			y = y(id(1:nid));
		else
			sxcell = num2cell(sx,1);
			y = spinterp(z, sxcell{:});
		end
	else
		% no points for this index set; continue with next 
		skip = 1;
		break;
	end
	
	if ~skip
		if isminimize
			[ymintemp, id] = min(y);
			if ymintemp < ymin
				ymin = ymintemp;
				xoptmin = sx(id,:);
			end
		end
		if ismaximize
			[ymaxtemp, id] = max(y);
			if ymaxtemp > ymax
				ymax = ymaxtemp;
				xoptmax = sx(id,:);
			end
		end
	end
end

xopt = [xoptmin' xoptmax'];
fval = [ymin ymax];
