function solver_output = odefinalize(obj,solver, sol,...
                                     outfun, outargs,...
                                     printstats, statvect,...
                                     nout,...
                                     haveeventfun, teout, yeout, ieout,...
                                     interp_data)
%ODEFINALIZE Helper function called by ODE solvers at the end of integration.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S, 
%            ODE23T, ODE23TB, ODE45, DDE23, DDESD.

%   Jacek Kierzenka
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2010/08/23 23:09:25 $

if ~isempty(outfun)
  feval(outfun,[],[],'done',outargs{:});
end

% Return more stats for implicit solvers: ODE15i, ODE15s, ODE23s, ODE23t, ODE23tb
fullstats = (length(statvect) > 3);  % faster than 'switch' or 'ismember'

stats = struct('nsteps',statvect(1),'nfailed',statvect(2),'nfevals',statvect(3)); 
if fullstats
  stats.npds     = statvect(4);
  stats.ndecomps = statvect(5);
  stats.nsolves  = statvect(6);  
else 
  statvect(4:6) = 0;   % Backwards compatibility
end  

if printstats
  fprintf('%g successful steps\n', stats.nsteps);
  fprintf('%g failed attempts\n', stats.nfailed);
  fprintf('%g function evaluations\n', stats.nfevals);
  if fullstats
    fprintf('%g partial derivatives\n', stats.npds);
    fprintf('%g LU decompositions\n', stats.ndecomps);
    fprintf('%g solutions of linear systems\n', stats.nsolves);
  end
end

solver_output = {};

if (nout > 0) % produce output
  if isempty(sol) % output [t,y,...]
    if nout < length(obj.t)
        obj.t=obj.t(1:nout);
        obj.y=obj.y(:,1:nout);
    end
    %solver_output{1} = obj.t(1:nout);
    %solver_output{2} = obj.y(:,1:nout);
    if haveeventfun
      solver_output{3} = teout;
      solver_output{4} = yeout;
      solver_output{5} = ieout;
    end
    solver_output{end+1} = statvect(:);  % Column vector
  else % output sol  
    % Add remaining fields
    sol.x = tout(1:nout);
    sol.y = yout(:,1:nout);
    if haveeventfun
      sol.xe = teout;
      sol.ye = yeout;
      sol.ie = ieout;
    end
    sol.stats = stats;
    switch solver
     case {'dde23','ddesd'}
      [history,ypout] = deal(interp_data{:});
      sol.yp = ypout(:,1:nout);
      if isstruct(history)
        sol.x = [history.x sol.x];
        sol.y = [history.y sol.y];
        sol.yp = [history.yp sol.yp];
        if isfield(history,'xe')
          if isfield(sol,'xe')
            sol.xe = [history.xe sol.xe];
            sol.ye = [history.ye sol.ye];
            sol.ie = [history.ie sol.ie];
          else
            sol.xe = history.xe;
            sol.ye = history.ye;
            sol.ie = history.ie;
          end
        end
      end
     case 'ode45'
      [f3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.f3d = f3d(:,:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode15s'      
      [kvec,dif3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.kvec = kvec(1:nout);
      maxkvec = max(sol.idata.kvec);
      sol.idata.dif3d = dif3d(:,1:maxkvec+2,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode113'
      [klastvec,phi3d,psi2d,idxNonNegative] = deal(interp_data{:});
      sol.idata.klastvec = klastvec(1:nout);
      kmax = max(sol.idata.klastvec);
      sol.idata.phi3d = phi3d(:,1:kmax+1,1:nout);
      sol.idata.psi2d = psi2d(1:kmax,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode23'
      [f3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.f3d = f3d(:,:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode23s'
      [k1data,k2data] = deal(interp_data{:});
      sol.idata.k1 = k1data(:,1:nout);
      sol.idata.k2 = k2data(:,1:nout);
     case 'ode23t'
      [zdata,znewdata,idxNonNegative] = deal(interp_data{:});
      sol.idata.z = zdata(:,1:nout);
      sol.idata.znew = znewdata(:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode23tb'
      [t2data,y2data,idxNonNegative] = deal(interp_data{:});
      sol.idata.t2 = t2data(1:nout);
      sol.idata.y2 = y2data(:,1:nout);           
      sol.idata.idxNonNegative = idxNonNegative;
     case 'ode15i'      
      [kvec,ypfinal] = deal(interp_data{:});
      sol.idata.kvec = kvec(1:nout);
      sol.extdata.ypfinal = ypfinal;
     otherwise
      error(message('MATLAB:odefinalize:UnrecognizedSolver', solver));
    end  
    if strcmp(solver,'dde23') || strcmp(solver,'ddesd')
      solver_output = sol;
    else  
      solver_output{1} = sol;
    end  
  end
end    
