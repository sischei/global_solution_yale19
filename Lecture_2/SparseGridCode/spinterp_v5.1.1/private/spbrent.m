function [xopt, fval, exitflag, output] = spbrent(z,pvec,xdir,xvec,fxvec,options)
% SPBRENT  Brent's line search algorithm variant for sparse grids
%    [XOPT, FVAL] = SPBRENT(Z,P,XDIR,X,FX)  Finds a local minimizer
%    along the line f(y) = P + XDIR*y for the bracket X = [a,b,c], 
%    FX = [f(a),f(b),f(c)], with a < b < c and f(a) > f(b) < f(c).
%    
%    [XOPT, FVAL] = SPBRENT(Z,P,XDIR,X,FX,OPTIONS)  Allows to
%    specify optimization options set with spoptimset.
%
%   [X,FVAL,EXITFLAG] = SPBRENT(...)  returns an EXITFLAG that 
%   describes the exit condition of SPBRENT. Possible values of 
%   EXITFLAG and the corresponding exit conditions are
%
%    1  SPBRENT converged to a solution X.
%    0  Maximum number of function evaluations or iterations
%       reached.
%
%
%    See also SPINTERPLINE, SPOPTIMSET
	
% Author : Andreas Klimke
% Version: 1.0
% Date   : September 19, 2006

% Change log:
% V1.0   : September 19, 2006
%          Initial version

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

t0 = clock;

if nargin < 6, options = []; end 

dispopt = spoptimget(options, 'Display', 'off');
[isdispiter, iterstr] = initoptidisp(dispopt);
tolx = spoptimget(options, 'TolX', 1e-4);
% Default is ignore TolFun
tolfun = spoptimget(options, 'TolFun', 0); 
maxiter = spoptimget(options, 'MaxIter', 100);
if isfield(z,'selectOutput')
	numout = z.selectOutout;
else
	numout = 1;
end	

% Absolute tolerance to avoid computing values close to zero beyond
% floating point accuracy
abstol = (z.fevalRange(numout,2) - z.fevalRange(numout,1)).*100*eps;

goldsec = 1 - (3 - sqrt(5)) / 2;
nfevals = 0;

% Initialize best points with bracketing triplet
if xvec(1) < xvec(3) 
	a = xvec(1);
	b = xvec(3);
	if fxvec(1) < fxvec(3)
		v = a; 
		fv = fxvec(1);
		w = b; 
		fw = fxvec(3);
	else
		v = b;
		fv = fxvec(3);
		w = a;
		fw = fxvec(1);
	end
else
	a = xvec(3);
	b = xvec(1);
	if fxvec(1) < fxvec(3)
		v = b;
		fv = fxvec(3);
		w = a;
		fw = fxvec(1);
	else
		v = a;
		fv = fxvec(1);
		w = b;
		fw = fxvec(3);
	end
end
x = xvec(2);
fx = fxvec(2);
if isdispiter
	disp(sprintf(iterstr, 0, 0, 0, fw, 'start point'));
	disp(sprintf(iterstr, 0, 0, 0, fv, 'start point'));
	disp(sprintf(iterstr, 0, 0, 0, fx, 'start point'));
end

% Initialize d
xm = 0.5*(a+b);
if x >= xm,	e = a - x; else e = b - x; end
d = e;
e = 0;

exitflag = 0;

% Main interation loop
for k = 1:maxiter
  tol1 = sqrt(eps) * abs(x) + tolx/3;
	tol2 = 2 * tol1;
	xm = 0.5 * (a+b);
	if 2.0*abs(fw-fx) <= max(tolfun * (abs(fw)+abs(fx)),abstol) || ...
			abs(x-xm) <= tol2 - 0.5*(b-a)
	  exitflag = 1;
		break;
	end
	if abs(e) > tol1
		r = (x - w) * (fx - fv);
		q = (x - v) * (fx - fw);
		p = (x - v) * q - (x - w)*r;
		q = 2 * (q - r);
		if q > 0
		  p = -p;
		end
		q = abs(q);
		etemp = e;
		e = d;
		if abs(p) >= abs(0.5 * q * etemp) || p <= q*(a - x) || p >= q*(b-x)
		  if x >= xm
			  e = a - x;
			else
			  e = b - x;
			end
		  d = goldsec * e;
  		how = 'golden';
		else
		  d = p / q;
			u = x + d;
			if (u - a < tol2) || (b - u < tol2)
			  d = tol1 * (sign(xm - x) + ((xm - x) == 0));
			end
			how = 'parabolic fit';
		end
	else
	  if x >= xm
		  e = a - x;
		else
		  e = b - x;
		end
	  d = goldsec * e;
		how = 'golden';
  end
	if abs(d) >= tol1
	  u = x + d;
	else
	  u = x + (sign(d) + (d == 0)) * tol1;
	end
	fu = spinterpline(z, pvec, xdir, u);
  nfevals = nfevals + 1;
  if isdispiter, disp(sprintf(iterstr, nfevals, nfevals, 0, fu, how)); end			
	if fu <= fx
	  if u >= x
		  a = x;
		else
		  b = x;
		end
		v = w; w = x; x = u;
		fv = fw; fw = fx; fx = fu;
	else
	  if u < x
		  a = u;
		else
		  b = u;
	  end
		if fu <= fw || w == x
		  v = w;
			w = u;
			fv = fw;
			fw = fu;
		elseif fu <= fv || v == x || v == w
		  v = u;
			fv = fu;
		end
	end
end
xopt = x;
fval = fx;

% Return stats
if nargout == 4
  output.nFEvals = nfevals;
	output.time = etime(clock, t0);
end
