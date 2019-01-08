function w = spquadwgp(levelseq,w1d,id)
% SPQUADWGP   Compute quadrature weights, Gauss-Patterson grid 
%    W = SPQUADWGP(LEVELSEQ,W1D,ID)  Computes the quadrature 
%    weights for the given sequence of index sets LEVELSEQ. One
%    row of column vectow W represents one weight. 
%    W1D are the integration weights in 1D, ID are the indices
%    where the weights of an according level start. 
%    (Internal function)
	
% Author : Andreas Klimke
% Version: 1.0
% Date   : November 24, 2007

% Change log:
% V1.1   : November 24, 2007
%          Initial version

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------


% Get the number of levels
nlevels = uint32(size(levelseq,1));

% Get the dimension
d = uint16(size(levelseq,2));

% Init weights vector
w = zeros(spgetnpointsnb(levelseq),1);
	
% index contains the index of the resulting array containing all
% subdomains of the level.
index = 1;

wid = zeros(d,1);
widstart = zeros(d,1);
widend = zeros(d,1);

for kl = 1:nlevels
	npoints = 1;
	for k = 1:d
		lval = double(levelseq(kl,k));
		np = 2^lval;
		npoints = npoints * np;
		widstart(k) = id(lval+1);
		wid(k) = widstart(k);
		widend(k) = wid(k) + np - 1;
	end
	
	wid1 = wid(1);
	for k = index:index+npoints-1
		
		% Compute multi-dimensional weight
		wval = w1d(wid1);
		for l = 2:d
			wval = wval * w1d(wid(l));
		end
		w(k) = wval;
		
		% Compute next set of indices into w1d
		if wid1 < widend(1)
			wid1 = wid1 + 1;
		else
			wid1 = widstart(1);
			for l = 2:d
				if wid(l) < widend(l)
					wid(l) = wid(l) + 1;
					break;
				else
					wid(l) = widstart(l);
				end
			end
		end
	end
		
	index = index + npoints;
end
