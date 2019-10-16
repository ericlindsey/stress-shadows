function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, args, odeFcn, ...
    options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
    dataType, verbose] =   ...
    odearguments(obj,FcnHandlesUsed, solver, ode, tspan, y0, options, extras)
%ODEARGUMENTS  Helper function that processes arguments for all ODE solvers.
%
%   See also ODE113, ODE15I, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Mike Karr, Jacek Kierzenka
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.12.4.15 $  $Date: 2011/05/17 02:23:32 $

if strcmp(solver,'ode15i')
    FcnHandlesUsed = true;   % no MATLAB v. 5 legacy for ODE15I
end

if FcnHandlesUsed  % function handles used
    if isempty(tspan) || isempty(y0)
        error(message('MATLAB:odearguments:TspanOrY0NotSupplied', solver));
    end
    if length(tspan) < 2
        error(message('MATLAB:odearguments:SizeTspan', solver));
    end
    htspan = abs(tspan(2) - tspan(1));
    tspan = tspan(:);
    ntspan = length(tspan);
    t0 = tspan(1);
    next = 2;       % next entry in tspan
    tfinal = tspan(end);
    args = extras;                 % use f(t,y,p1,p2...)
    
else  % ode-file used   (ignored when solver == ODE15I)
    % Get default tspan and y0 from the function if none are specified.
    if isempty(tspan) || isempty(y0)
        if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 )
            error(message('MATLAB:odearguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));
        end
        [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
        if isempty(tspan)
            tspan = def_tspan;
        end
        if isempty(y0)
            y0 = def_y0;
        end
        options = odeset(def_options,options);
    end
    tspan = tspan(:);
    ntspan = length(tspan);
    if ntspan == 1    % Integrate from 0 to tspan
        t0 = 0;
        next = 1;       % Next entry in tspan.
    else
        t0 = tspan(1);
        next = 2;       % next entry in tspan
    end
    htspan = abs(tspan(next) - t0);
    tfinal = tspan(end);
    
    % The input arguments of f determine the args to use to evaluate f.
    if (exist(ode)==2)
        if (nargin(ode) == 2)
            args = {};                   % f(t,y)
        else
            args = [{''} extras];        % f(t,y,'',p1,p2...)
        end
    else  % MEX-files, etc.
        try
            args = [{''} extras];        % try f(t,y,'',p1,p2...)
            feval(ode,tspan(1),y0(:),args{:});
        catch ME
            args = {};                   % use f(t,y) only
        end
    end
end

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if t0 == tfinal
    error(message('MATLAB:odearguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
    error(message('MATLAB:odearguments:TspanNotMonotonic'));
end

f0 = feval(ode,t0,y0,args{:});   % ODE15I sets args{1} to yp0.
[m,n] = size(f0);
if n > 1
    error(message('MATLAB:odearguments:FoMustReturnCol', funstring( ode )));
elseif m ~= neq
    error(message('MATLAB:odearguments:SizeIC', funstring( ode ), m, neq, funstring( ode )));
end

% Determine the dominant data type
classT0 = class(t0);
classY0 = class(y0);
classF0 = class(f0);
if strcmp(solver,'ode15i')
    classYP0 = class(args{1});  % ODE15I sets args{1} to yp0.
    dataType = superiorfloat(t0,y0,args{1},f0);
    
    if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
            strcmp(classF0,dataType) && strcmp(classYP0,dataType))
        input1 = '''t0'', ''y0'', ''yp0''';
        input2 = '''f(t0,y0,yp0)''';
        warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
    end
else
    dataType = superiorfloat(t0,y0,f0);
    
    if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
            strcmp(classF0,dataType))
        input1 = '''t0'', ''y0''';
        input2 = '''f(t0,y0)''';
        warning(message('MATLAB:odearguments:InconsistentDataType',input1,input2,solver));
    end
end

% Get the error control options, and set defaults.
verbose = unicycle.ode.odeget(options,'Verbose',true);
rtol = unicycle.ode.odeget(options,'RelTol',1e-3,'fast');
if (length(rtol) ~= 1) || (rtol <= 0)
    error(message('MATLAB:odearguments:RelTolNotPosScalar'));
end
if rtol < 100 * eps(dataType)
    rtol = 100 * eps(dataType);
    warning(message('MATLAB:odearguments:RelTolIncrease', sprintf( '%g', rtol )))
end
atol = unicycle.ode.odeget(options,'AbsTol',1e-6,'fast');
if any(atol <= 0)
    error(message('MATLAB:odearguments:AbsTolNotPos'));
end
normcontrol = strcmp(unicycle.ode.odeget(options,'NormControl','off','fast'),'on');
if normcontrol
    if length(atol) ~= 1
        error(message('MATLAB:odearguments:NonScalarAbsTol'));
    end
    normy = norm(y0);
else
    if (length(atol) ~= 1) && (length(atol) ~= neq)
        error(message('MATLAB:odearguments:SizeAbsTol', funstring( ode ), neq));
    end
    atol = atol(:);
    normy = [];
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = min(abs(tfinal-t0), abs(unicycle.ode.odeget(options,'MaxStep',0.1*(tfinal-t0),'fast')));
if hmax <= 0
    error(message('MATLAB:odearguments:MaxStepLEzero'));
end
htry = abs(unicycle.ode.odeget(options,'InitialStep',[],'fast'));
if ~isempty(htry) && (htry <= 0)
    error(message('MATLAB:odearguments:InitialStepLEzero'));
end

odeFcn = ode;

end