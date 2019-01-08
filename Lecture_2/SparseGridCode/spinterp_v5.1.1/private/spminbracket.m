function [xbrac, fxbrac, bflag, x, fx, fxgrad, nfevals] = ...
  spminbracket(z,x,fx,fxgrad,xdir,xbox,options,stepsize)
% Find a bracket enclosing a local minimizer (internal function) 

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------

if nargin < 7, options = []; end
if nargin < 8, stepsize = []; end

dispopt = spoptimget(options, 'Display', 'off');
[isdispiter, iterstr] = initoptidisp(dispopt);

maxfactor = 100;

absxdir = norm(xdir);
if absxdir == 0.0
	error('MATLAB:spinterp:invalidzero',['Search direction for ' ...
	  'spminbracket routine is an all-zero vector.']);
end

d = z.d;
if isempty(z.range)
	z.range = [zeros(d,1) ones(d,1)];
end
if isempty(xbox)
	xbox = z.range;
end

% Get the maximum sparse grid levels
if isfield(z, 'dimAdapt')
	n = z.maxLevel';
else
	n = z.maxLevel*ones(d,1,'uint8');
end
goldsec = 2 - (3 - sqrt(5)) / 2;

% Compute normalized direction; make sure we descend by computing
% sign of the directional derivative

dirsign = sign(-dot(fxgrad,xdir));
nxdir = xdir/absxdir * dirsign;

if isempty(stepsize)
	% Compute combined initial step size (with normalized direction 
	% vector)
	% Compute step size in each dimension
	dx = 1./2.^(double(n)) .* (z.range(:,2)-z.range(:,1));
	stepsize = abs(dot(dx,nxdir));
end

success = 0;
nfevals = 0;

% Convert to point on line
p = x;
xa = 0;
fxa = fx;
if isdispiter, disp(sprintf(iterstr, 0, 0, 0, fx, 'start point')); end

% Get next point; ensure it is inside search box
xb = stepsize;
[xb, bflag, fxb, fxbgrad] = evalCheckBorder(z,xb,p,nxdir,xbox);
niter = 1;
nfevals = nfevals + 1;

% The next point lies inside the search box
if bflag == 0
	if isdispiter, disp(sprintf(iterstr, 1, 1, 0, fxb, 'initial')); end
	% Function value is greater than initial one -> flip search direction
	if fxb > fxa
		xc = xa; xa = xb; xb = xc;
		fxc = fxa; fxa = fxb; fxb = fxc;
		stepsize = -stepsize;
	end

	stepsize = stepsize * goldsec;
	% Get third initial point
	xc = xb + stepsize;
	[xc, bflag, fxc, fxcgrad] = evalCheckBorder(z,xc,p,nxdir,xbox);
	nfevals = nfevals + 1;
	if bflag == 0
		if isdispiter, disp(sprintf(iterstr, niter, nfevals, 0, fxc, 'golden')); end
	elseif bflag == 1 && fxc <= fxb
		% The boundary point is the minimum
		xbrac = [];
		fxbrac = [];
		x = p + nxdir * xc;
		fx = fxc;
		fxgrad = fxcgrad;
		if isdispiter, disp(sprintf(iterstr, niter, nfevals, 1, fxc, 'boundary minimum')); end
		% return immediately
		return;
	else
	  bflag = 2; % could be 1 if fxb < fxc
		if isdispiter, disp(sprintf(iterstr, niter, nfevals, 1, fxc, 'boundary')); end
		success = 1;
	end
elseif bflag == 1 && fxb <= fxa
	% The boundary point is the minimum
	xbrac = [];
	fxbrac = [];
	x = p + nxdir * xb;
	fx = fxb;
	fxgrad = fxbgrad;
	if isdispiter, disp(sprintf(iterstr, niter, nfevals, 1, fxb, 'boundary minimum')); end
	% return immediately
	return;
else
  bflag = 2; % could be 1 if fxa < fxb
	if isdispiter, disp(sprintf(iterstr, niter, nfevals, 1, fxb, 'boundary')); end
	% The new point lies on the boundary, but it is not the minimum
	% Thus, get the third point. 
	xc = stepsize * (2 - goldsec);
	fxc = spinterpline(z,p,nxdir,xc);
	nfevals = nfevals + 1;
	% We are done, since we know there is an internal minimum
	% somewhere (but we have not bracketed it)
	if isdispiter, disp(sprintf(iterstr, niter, nfevals, 0, fxc, 'golden')); end
	success = 1;
end

if success ~= 1
	parafit = 1;
	while fxc < fxb % the third point is smaller than the second, this means
			            % we have not yet bracketed
		
		niter = niter + 1;
		if parafit == 1
			parafit = 1; % default is try parabolic fit in next step
			% Compute u by parabolic extrapolation; 
			r = (xb - xa) * (fxb - fxc);
			q = (xb - xc) * (fxb - fxa);
			u = xb - ((xb - xc) * q - (xa - xb) * r) / ...
					(2 * max(abs(q-r),eps) * (sign(q-r) + ((q-r) == 0)));
			ulim = xb + maxfactor * (xc - xb);
			if (xb - u) * (u - xc) > 0 
				% Parabolic u between b and c; try it 
				% No need to check for boundary since point is between
				% b and c
				fu = spinterpline(z,p,nxdir,u);
				nfevals = nfevals + 1;
				if isdispiter, disp(sprintf(iterstr, niter, nfevals, 0, fu, 'parabolic fit')); end
				if fu < fxc
					% Got a minimum between b and c
					xa = u;
					fxa = fu;
					break;
				elseif fu > fxb
					% Got a minimum between a and u
					xc = u;
					fxc = fu;
					break;
				else
					% Parabolic fit was no use, since fxb > fxu > fxc.
					% No new information, take golden section step
					u = xc + goldsec * (xc - xb);
					how = 'golden';
				end
			elseif (xc - u) * (u - ulim) > 0
				% Parabolic fit is between c and allowable limit
				% Do not allow parabolic fit if the current fit fails to 
				% produce a bracketing triplet
				how = 'parabolic fit';
				parafit = 0;
			elseif (u - ulim) * (ulim - xc) > 0
				% Allowable limit is exceeded; limit u
				how = 'parabolic limit';
				u = ulim;
			else
				% Reject parabolic u and do golden section
				how = 'golden';
				u = xc + goldsec * (xc - xb); 
			end
		else % parafit == 0
			% Do golden section
			% Next step, try again a parabolic fit.
			parafit = 1;
			how = 'golden';
			u = xc + goldsec * (xb - xa);
		end
		% If we get here, there is an u prepared to be evaluated
		[u, bflag, fu, fugrad] = evalCheckBorder(z,u,p,nxdir,xbox);
		
		nfevals = nfevals + 1;

		if bflag == 1
			if isdispiter, disp(sprintf(iterstr, niter, nfevals, 1, fu, 'boundary minimum')); end
			% The boundary point is the minimum
			xbrac = [];
			fxbrac = [];
			x = p + nxdir * u;
			fx = fu;
			fxgrad = fugrad;
			% return immediately
			return;
		else
			% Eliminate oldest point
			xa = xb;   xb = xc;   xc = u;
			fxa = fxb; fxb = fxc; fxc = fu;
			if bflag == 2
				if isdispiter, disp(sprintf(iterstr, niter, nfevals, 1, fu, 'boundary')); end
				% The new point lies on the boundary, but it is not the minimum
				% We are done, since we know there is an internal minimum
				% somewhere (but we have not bracketed it)
				break;
			else
				if isdispiter, disp(sprintf(iterstr, niter, nfevals, 0, fu, how)); end			
			end
		end
	end				
end

% Return the bracketing triplet in sorted order
[xbrac, id] = sort([xa xb xc]/absxdir*dirsign);  
fxbrac = [fxa, fxb, fxc];
fxbrac = fxbrac(id);
x = [];
fx = [];
fxgrad = [];

return;
	
% ------------------------------------------------------------------
function [x, bflag, fxb, fxbgrad] = evalCheckBorder(z,x,p,xdir,xbox)

xvec = p + x * xdir;
xsign = sign(x);

% Check if point outside of the sparse grid range
borderleft  = xvec - xbox(:,1) < 0; 
borderright = xbox(:,2) - xvec < 0;
if any(borderleft) || any(borderright)
	% Must intersect with the boundary of the range cube

	% Get the active boundary component; 
	maxr = 0;
	for k = 1:length(xvec)
		if xdir(k) == 0, continue; end
		if xsign * xdir(k) < 0
			if (p(k) == xbox(k,1)), 
			  r = inf;
			else
				r = (p(k)-xvec(k)) / (p(k)-xbox(k,1));
			end
		else
			if (p(k) == xbox(k,2))
			  r = inf; 
			else
				r = (xvec(k)-p(k)) / (xbox(k,2)-p(k));
			end
		end
		if r > maxr
			maxr = r;
			id = k;
		end
	end
		
	if xsign * xdir(id) < 0
		x = (xbox(id,1) - p(id)) / xdir(id);
	else
		x = (xbox(id,2) - p(id)) / xdir(id);
	end
	
	% Compute new point on the boundary of the box
	xvec = p + x * xdir;

	% We are already done, but we must check if boundary point is a
	% minimum. To check this, we compute the gradient, and see if
	% the directional derivative points toward the boundary for the
	% given search direction. 
	% Negative xsign means we are at the left bound w.r.t. the search
	% direction xdir, so the directional derivative should be
	% positive. Positive x means we are at the right bound w.r.t. the
	% search direction. Then, the directional derivative should
	% be negative. Thus, the sign of the directional derivative
	% must be the opposite of xsign.
	[fxb, fxbgrad] = spsurfun(xvec,z);
	if sign(dot(fxbgrad,xdir)) ~= sign(xsign) 
		% Minimum lies at the boundary. We can return this point and
		% its gradient. No line search for the minimum is required.
		bflag = 1;
	else
		% The sign of the directional derivative has changed; the 
		% minimum lies between the previous point and the boundary.
		bflag = 2;
	end
else
	fxb = spsurfun(xvec,z);
	fxbgrad = [];
	bflag = 0;	
end

