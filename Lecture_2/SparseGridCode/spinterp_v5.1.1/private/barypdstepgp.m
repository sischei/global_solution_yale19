function ip = barypdstepgp(z, allnx, dims, x, y, w)
% BARYPDSTEPGP  Step of Barycentric d-dim. Lagrange interpolation 
%    Gauss-Patterson grid. Includes additional weight vector w
%    containing barycentric weights
%    (internal function)

% Author : Andreas Klimke
% Version: 1.0
% Date   : December 4, 2007

% Change log:
% V1.0   : Initial release.

% ------------------------------------------------------------
% Sparse Grid Interpolation Toolbox
% Copyright (c) 2006 W. Andreas Klimke, Universitaet Stuttgart 
% Copyright (c) 2007-2008 W. A. Klimke. All Rights Reserved.
% See LICENSE.txt for license. 
% email: klimkeas@ians.uni-stuttgart.de
% web  : http://www.ians.uni-stuttgart.de/spinterp
% ------------------------------------------------------------
	
d = uint8(length(allnx));

nt = uint32(size(y,1));
idaddfix = zeros(d,1,'uint32');
idadd = zeros(d,1,'uint32');
aidstart = ones(d,1,'uint32');
aidmax = zeros(d,1,'uint32');
aidvec = zeros(d,1,'uint32');

ip = zeros(nt,1);
c1 = zeros(d,1);
c2 = zeros(d,1);
id = zeros(d,1,'uint32');
nx = zeros(d,1,'uint32');
idvec = zeros(d,1,'uint32');
c11 = 0;
c12 = 0;

aid = uint32(0);
for k1 = 1:d
  % modify nx such that it only contains the non-zero supports
  nx(k1) = (allnx(k1)+1)/2;
end

a = zeros(sum(nx),1);

idaddfix(1) = 1;
for k2 = 2:d
  idaddfix(k2) = idaddfix(k2-1)*nx(k2-1);
  aidstart(k2) = aidstart(k2-1) + nx(k2-1);
end

for k3 = 1:d
  aidmax(k3) = (nx(k3)+aidstart(k3))-1;
end

for m = 1:nt

  xid = uint32(0);
  aid = uint32(0);
  iszeronode = 0;
  
  for l = 1:d
    t = y(m,dims(l));
    frac = 0;
    s = 0;
    id(l) = 0;
    actnx = allnx(l);
    actna = nx(l);
        
    if actnx == 3
      if t == x(xid+1)
        id(l) = 1;
        c2(l) = 1;
        aid = aid + 2;
        xid = xid + 3;
        continue;
      elseif t == x(xid+2)
        ip(m) = 0;
				iszeronode = 1;
        break;
      elseif t == x(xid+3)
        id(l) = 2;
        c2(l) = 1;
        aid = aid + 2;
        xid = xid + 3;
        continue;
      end
      
      % do first element
      frac = w(xid+1) / (t-x(xid+1));
      a(aid+1) = frac;
      s = s + frac;
    
      % do second element
      s = s + w(xid+2) / (t-x(xid+2));
      
      % do third element
      frac = w(xid+3) / (t-x(xid+3));
      a(aid+2) = frac;
      s = s + frac;  
    else
      if t >= 0 && t <= 1
        k = uint32(xid+actnx/2);
        lok = uint32(xid+1);
        upk = uint32(xid+actnx);
				xk = x(k);
        while lok + 1 < upk
          if xk > t 
            upk = k;
          elseif xk <= t
            lok = k;
          end
					k = (lok+upk)/2;
          xk = x(k);
        end
        if t == x(lok)
          k = lok;
        elseif t == x(upk)
          k = upk;
        else
          k = uint32(0);
        end
        if k > 0
          if mod(k-xid,2) == 0
            iszeronode = 1;
            ip(m) = 0;
            break;
          else
            id(l) = (k-xid+1)/2;
            c2(l) = 1;
            aid = aid + actna;
            xid = xid + actnx;
            continue;
          end
        end
      end
    
      % do every element
      kxi = xid + 1;
      for k = aid+1:aid+actna-1
        frac = w(kxi)/(t-x(kxi));
        a(k) = frac;
        s = s + frac + w(kxi+1)/(t-x(kxi+1));
        kxi = kxi + 2;
      end
      
      % do last element
      frac = w(kxi)/(t-x(kxi));
      a(aid+actna) = frac;
      s = s + frac;
      
    end
    c2(l) = s;
    aid = aid+actna;
    xid = xid+actnx;
  end
  if iszeronode, continue; end
  
  niter = uint32(1);
  nid = uint8(0);
  for l = 1:d
    if id(l) == 0
      niter = niter * nx(l);
      nid = nid + 1;
      idvec(nid) = l;
    else
      % Modify id such that the first index is 0, not 1. But
      % otherwise, leave id as it is; it remains constant
      % throughout the interpolation process.
      id(l) = id(l) - 1;
    end
    % disp([sprintf('id[%d] = %d', l, id(l))]);
  end
  if nid == 0
    % Point to interpolate lies on a grid point - set result to
    % this value.
    zid = uint32(1);
    for l = 1:d
      zid = zid + id(l)*idaddfix(l);
    end
    ip(m) = z(zid);
    % disp(sprintf('ip(%d) = %f',m,ip(m)));
  else
    zid = uint32(1);
    for l = 1:d
      c1(l) = 0;
      zid = zid + id(l)*idaddfix(l);
      aidvec(l) = aidstart(l);
    end
    
    for l = 1:nid
      % modify idadd such that it correctly adds only the
      % difference from the previous z index. This must take into
      % account which dimensions are skipped in case the point to
      % interpolate lies on some grid axes.
      
      % get the next active index
      idl = idvec(l);
      if l > 1
        idadd(idl) = idaddfix(idl)-idaddfix(idvec(l-1)+1);
        % idadd(idl) becomes zero of the previous index was not
        % skipped. 
      else
        idadd(idl) = idaddfix(idl);
      end
    end

    c11 = 0;
    c12 = 0;

    id1 = idvec(1);
    addid1 = idadd(id1);
    aid1 = aidvec(id1);
    aidmax1 = aidmax(id1);
    
    if nid == 1
      for k2 = aid1:aidmax1
        c11 = c11 + z(zid)*a(k2);
        zid = zid + addid1;
      end
      c1d = c11;
    else
      niter = niter/nx(id1);
      id2 = idvec(2);
      aid2 = aidvec(id2);
      aidmax2 = aidmax(id2);
      addid2 = idadd(id2);
      if nid == 2
        k = uint32(0);
        while k < niter
          for k2 = aid1:aidmax1
            c11 = c11 + z(zid)*a(k2);
            zid = zid + addid1;
          end
          zid = zid + addid2;
          c12 = c12 + a(aid2) * c11;
          c11 = 0;
          aid2 = aid2 + 1;
          k = k + 1;
        end
        c1d = c12;
      else
        % Interpolation in three or more dimensions
        k = uint32(0);
        while k < niter
          for k2 = aid1:aidmax1
            c11 = c11 + z(zid)*a(k2);
            zid = zid + addid1;
          end
          zid = zid + addid2;
          c12 = c12 + a(aid2) * c11;
          c11 = 0;
          if aid2 == aidmax2
            zid = zid + idadd(idvec(3));
            c1(3) = c1(3) + a(aidvec(idvec(3))) * c12;
            c12 = 0;
            aid2 = aidstart(id2);
            for l = 3:nid
              idl = idvec(l);
              aidl = aidvec(idl);
              if aidl == aidmax(idl)
                if l < nid
                  zid = zid + idadd(idvec(l+1));
                  c1(l+1) = c1(l+1) + a(aidvec(idvec(l+1))) * c1(l);
                  c1(l) = 0;
                  aidvec(idl) = aidstart(idl);
                end
              else
                aidvec(idl) = aidl + 1;
                break;
              end
            end
          else
            aid2 = aid2 + 1;
          end
          k = k + 1;
        end
        c1d = c1(nid);
      end      
    end
    
    c2d = 1;
    for l = 1:d
      c2d = c2d * c2(l);
    end
    ip(m) = c1d/c2d;
  end
end
