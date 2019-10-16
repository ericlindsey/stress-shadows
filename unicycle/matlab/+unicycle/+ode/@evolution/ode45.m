function varargout = ode45(obj,tspan,y0,options,varargin)
% ODE45 solve ordinary differential equation defined by the differential
% operator evolution.ode(t,y), using a 4-order accurate method and automatic
% time steps.
%
%   [TOUT,YOUT] = UNICYCLE.ODE.EVOLUTION/ODE45(TSPAN,Y0,OPTIONS)
%
% ODE45 solve the system of differential equations y' = f(t,y) from time T0
% to TFINAL with initial conditions Y0.
%
% INPUT:
%
%   TSPAN - period of integration [T0 TFINAL]
%
% OUTPUT:
%
%   TOUT  - vector of epochs when solutions are available
%   YOUT  - matrix of solutions (one line per time step)
%
%   Each row in the solution array YOUT corresponds to a time returned in
%   the column vector TOUT. To obtain solutions at specific times
%   T0,T1,...,TFINAL (all increasing or all decreasing), use
%
%     TSPAN = [T0 T1 ... TFINAL].
%
% To alter the integration properties, use:
%
%   [TOUT,YOUT] = EVOLUTION/ODE45(TSPAN,Y0,OPTIONS)
%
% where OPTIONS created with the ODESET function. See ODESET for details.
%
% To retrieve solutions at particular states (events), use:
%
%   [TOUT,YOUT,TE,YE,IE] = EVOLUTION/ODE45(ODEFUN,TSPAN,Y0,OPTIONS)
%
% with the 'Events' property in OPTIONS set to a function handle EVENTS,
% which finds where functions of (T,Y), called event functions, are zero.
% For each function you specify whether the integration is to terminate at
% a zero and whether the direction of the zero crossing matters. These are
% the three column vectors returned by EVENTS:
%
%   [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y).
%
% For the I-th event function:
%
% VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration
% is to terminate at a zero of this event function and 0 otherwise.
% DIRECTION(I)=0 if all zeros are to be computed (the default), +1 if only
% zeros where the event function is increasing, and -1 if only zeros where
% the event function is decreasing. Output TE is a column vector of times
% at which events occur. Rows of YE are the corresponding solutions, and
% indices in vector IE specify which event occurred.
%
%   SOL = EVOLUTION/ODE45(ODEFUN,[T0 TFINAL],Y0...)
%
% returns a structure that can be used with DEVAL to evaluate the solution
% or its first derivative at any point between T0 and TFINAL. The steps
% chosen by ODE45 are returned in a row vector SOL.x.  For each I, the
% column SOL.y(:,I) contains the solution at SOL.x(I). If events were
% detected, SOL.xe is a row vector of points at which events occurred.
% Columns of SOL.ye are the corresponding solutions, and indices in vector
% SOL.ie specify which event occurred.
%
% ODE45 is an implementation of the explicit Runge-Kutta (4,5) pair of
% Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.

% Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
% M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

% Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
% Copyright 1984-2011 The MathWorks, Inc.
% $Revision: 5.74.4.13 $  $Date: 2011/04/16 06:38:58 $

import unicycle.ode.*
import unicycle.utils.*

solver_name = 'ode45';

% Check inputs
if nargin < 4
    options = [];
    if nargin < 3
        y0 = [];
        if nargin < 2
            tspan = [];
            if nargin < 1
                error(message('MATLAB:ode45:NotEnoughInputs'));
            end
        end
    end
end

ode=@obj.ode;

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0;

% Output
FcnHandlesUsed  = isa(ode,'function_handle');
output_sol = (FcnHandlesUsed && (nargout==1));      % sol = odeXX(...)
%output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)
output_ty = true;
% There might be no output requested...

% Coseismic events
evtIndex = 1;

sol = []; f3d = [];
if output_sol
    sol.solver = solver_name;
    sol.extdata.odefun = ode;
    sol.extdata.options = options;
    sol.extdata.varargin = varargin;
end

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
    options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
    dataType, verbose] = ...
    obj.odearguments(FcnHandlesUsed,solver_name,ode,tspan,y0,options,varargin);
nfevals = nfevals + 1;



% Handle the output
%if nargout > 0
outputFcn = unicycle.ode.odeget(options,'OutputFcn',[],'fast');
%else
%    outputFcn = unicycle.ode.odeget(options,'OutputFcn',@odeplot,'fast');
%end
outputArgs = {};
if isempty(outputFcn)
    haveOutputFcn = false;
else
    haveOutputFcn = true;
    outputs = unicycle.ode.odeget(options,'OutputSel',1:neq,'fast');
    if isa(outputFcn,'function_handle')
        % With MATLAB 6 syntax pass additional input arguments to outputFcn.
        outputArgs = varargin;
    end
end

refine = max(1,unicycle.ode.odeget(options,'Refine',4,'fast'));
if isfield(options,'OutputSelection')
    ioselection = options.OutputSelection;
else
    ioselection = [];
end
if isa(ioselection ,'function_handle')
    outputAt = 'Subset';
else
    if ntspan > 2
        outputAt = 'RequestedPoints';         % output only at tspan points
    elseif refine <= 1
        outputAt = 'SolverSteps';             % computed points, no refinement
    else
        outputAt = 'RefinedSteps';            % computed points, with refinement
        S = (1:refine-1) / refine;
    end
end

printstats = strcmp(unicycle.ode.odeget(options,'Stats','off','fast'),'on');
oDir = unicycle.ode.odeget(options,'oDir','','fast');

% Handle the event function
[haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
    obj.odeevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);

% Non-negative solution components
idxNonNegative = unicycle.ode.odeget(options,'NonNegative',[],'fast');
nonNegative = ~isempty(idxNonNegative);
if nonNegative  % modify the derivative function
    [odeFcn,thresholdNonNegative] = obj.odenonnegative(odeFcn,y0,threshold,idxNonNegative);
    f0 = feval(odeFcn,t0,y0,odeArgs{:});
    nfevals = nfevals + 1;
end

if verbose
    textprogressbar('# solving ode             : ');
end

t = t0;
y = y0;
yold= y0;

% Allocate memory if we're generating output.
nout = 0;
obj.t = []; obj.y = [];

chunkPoints = max(100,50*refine);
% initiate and allocate memory for observation points
if isobject(obj.flt)
    for k=1:length(obj.flt.observationPoints)
        obj.flt.observationPoints{k}.t=zeros(1,chunkPoints,dataType);
        obj.flt.observationPoints{k}.y=zeros(obj.flt.dgf,chunkPoints,dataType);
        obj.flt.observationPoints{k}.N=0;
    end
end

if isobject(obj.shz)
    for k=1:length(obj.shz.observationPoints)
        obj.shz.observationPoints{k}.t=zeros(1,chunkPoints,dataType);
        obj.shz.observationPoints{k}.y=zeros(obj.shz.dgf,chunkPoints,dataType);
        obj.shz.observationPoints{k}.N=0;
    end
end

% initiate and allocate memory for vMax and tMax
nMax=1;

if isobject(obj.flt)
    obj.flt.tMax=zeros(1,chunkPoints,dataType);
    obj.flt.vMax=zeros(1,chunkPoints,dataType);
    obj.flt.tMax(1)=t;
end

if isobject(obj.shz)
    obj.shz.tMax=zeros(1,chunkPoints,dataType);
    obj.shz.eMax=zeros(1,chunkPoints,dataType);
    obj.shz.tMax(1)=t;
end

if ~isempty(oDir)
    
    timeFID=fopen([oDir '/time.dat'],'wt');
    fprintf(timeFID,'# t(yr) dt(yr)\n');
    fprintf(timeFID,'%20.16e %20.16e\n',[0 0]);
    
    if isobject(obj.flt)
        vMaxFID=fopen([oDir '/vmax.dat'],'wt');
        fprintf(vMaxFID,'# time (yr), max speed (m/yr)\n');
        fprintf(vMaxFID,'%20.16e %20.16e\n',[0 0]);
    
        fltOpts=cell(length(obj.flt.observationPoints),1);
        for k=1:length(obj.flt.observationPoints)
            fltOpts{k}=fopen(sprintf('%s/flt-opts-%03d.dat',oDir,k),'wt');
            fprintf(timeFID,'# degrees of freedom\n');
        end
    end
    
    if isobject(obj.shz)
        eMaxFID=fopen([oDir '/emax.dat'],'wt');
        fprintf(eMaxFID,'# time (yr), max strain rate (1/yr)\n');
        fprintf(eMaxFID,'%20.16e %24.20e\n',[0 0]);
        
        shzOpts=cell(length(obj.shz.observationPoints),1);
        for k=1:length(obj.shz.observationPoints)
            shzOpts{k}=fopen(sprintf('%s/shz-opts-%03d.dat',oDir,k),'wt');
            fprintf(timeFID,'# degrees of freedom\n');
        end
    end
end

if isobject(obj.flt)
    % initiate fault observation points
    for k=1:length(obj.flt.observationPoints)
        obj.flt.observationPoints{k}.t(1)=t;
        obj.flt.observationPoints{k}.y(:,1)=y(1+(obj.flt.observationPoints{k}.index-1)*obj.flt.dgf:obj.flt.observationPoints{k}.index*obj.flt.dgf,1);
        obj.flt.observationPoints{k}.N=1;
    end
end

if isobject(obj.shz)
    % initiate shear zones observation points
    for k=1:length(obj.shz.observationPoints)
        obj.shz.observationPoints{k}.t(1)=t;
        obj.shz.observationPoints{k}.y(:,1)=y(1+(obj.shz.observationPoints{k}.index-1)*obj.shz.dgf:obj.shz.observationPoints{k}.index*obj.shz.dgf,1);
        obj.shz.observationPoints{k}.N=1;
    end
end

if output_sol
    chunk = min(max(100,50*refine), refine+floor((2^11)/neq));
    obj.t = zeros(1,chunk,dataType);
    obj.y = zeros(neq,chunk,dataType);
    f3d  = zeros(neq,7,chunk,dataType);
else
    if ntspan > 2                         % output only at tspan points
        obj.t = zeros(1,ntspan,dataType);
        obj.y = zeros(neq,ntspan,dataType);
    else                                  % alloc in chunks
        chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
        obj.t = zeros(1,chunk,dataType);
        obj.y = zeros(neq,chunk,dataType);
    end
end
nout = 1;
obj.t(nout) = t;
obj.y(:,nout) = y;

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1];
B = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
f = zeros(neq,7,dataType);

hmin = 16*eps(t);
if isempty(htry)
    % Compute an initial step size h using y'(t).
    absh = min(hmax, htspan);
    if normcontrol
        rh = (norm(f0) / max(normy,threshold)) / (0.8 * rtol^pow);
    else
        rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
    end
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
else
    absh = min(hmax, max(hmin, htry));
end
f(:,1) = f0;

% Initialize the output function.
if haveOutputFcn
    feval(outputFcn,[t tfinal],y(outputs),'init',outputArgs{:});
end

% THE MAIN LOOP

done = false;
while ~done
    
    % By default, hmin is a small number such that t+hmin is only slightly
    % different than t.  It might be 0 if t is 0.
    hmin = 16*eps(t);
    absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
    h = tdir * absh;
    
    % Stretch the step if within 10% of tfinal-t.
    if 1.1*absh >= abs(tfinal - t)
        h = tfinal - t;
        absh = abs(h);
        done = true;
    end
    
    if evtIndex<=length(obj.evt)
        if t ~= obj.evt{evtIndex}.src.t0
            if 1.1*absh >= abs(obj.evt{evtIndex}.src.t0 - t)
                h = obj.evt{evtIndex}.src.t0 - t;
                absh = abs(h);
            end
        end
    end
    
    
    
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true
        hA = h * A;
        hB = h * B;
        
        f(:,2) = feval(odeFcn,t+hA(1),y+f*hB(:,1),odeArgs{:});
        f(:,3) = feval(odeFcn,t+hA(2),y+f*hB(:,2),odeArgs{:});
        f(:,4) = feval(odeFcn,t+hA(3),y+f*hB(:,3),odeArgs{:});
        f(:,5) = feval(odeFcn,t+hA(4),y+f*hB(:,4),odeArgs{:});
        f(:,6) = feval(odeFcn,t+hA(5),y+f*hB(:,5),odeArgs{:});
        tnew = t + hA(6);
        if done
            tnew = tfinal;   % Hit end point exactly.
        end
        h = tnew - t;      % Purify h.
        
        ynew = y + f*hB(:,6);
        f(:,7) = feval(odeFcn,tnew,ynew,odeArgs{:});
        nfevals = nfevals + 6;
        
        % Estimate the error.
        NNrejectStep = false;
        if normcontrol
            normynew = norm(ynew);
            errwt = max(max(normy,normynew),threshold);
            err = absh * (norm(f * E) / errwt);
            if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
                if errNN > rtol
                    err = errNN;
                    NNrejectStep = true;
                end
            end
        else
            err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
            if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);
                if errNN > rtol
                    err = errNN;
                    NNrejectStep = true;
                end
            end
        end
        
        % Accept the solution only if the weighted error is no more than the
        % tolerance rtol.  Estimate an h that will yield an error of rtol on
        % the next step or the next try at taking this step, as the case may be,
        % and use 0.8 of this value to avoid failures.
        if err > rtol                       % Failed step
            nfailed = nfailed + 1;
            if absh <= hmin
                warning(message('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
                solver_output = obj.odefinalize(solver_name, sol,...
                    outputFcn, outputArgs,...
                    printstats, [nsteps, nfailed, nfevals],...
                    nout, obj.t, obj.y,...
                    haveEventFcn, teout, yeout, ieout,...
                    {f3d,idxNonNegative});
                if nargout > 0
                    varargout = solver_output;
                end
                return;
            end
            
            if nofailed
                nofailed = false;
                if NNrejectStep
                    absh = max(hmin, 0.5*absh);
                else
                    absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
                end
            else
                absh = max(hmin, 0.5 * absh);
            end
            h = tdir * absh;
            done = false;
            
        else                                % Successful step
            NNreset_f7 = false;
            if nonNegative && any(ynew(idxNonNegative)<0)
                ynew(idxNonNegative) = max(ynew(idxNonNegative),0);
                if normcontrol
                    normynew = norm(ynew);
                end
                NNreset_f7 = true;
            end
            
            break;
            
        end
    end
    nsteps = nsteps + 1;
    
    if evtIndex<=length(obj.evt)
        if t==obj.evt{evtIndex}.src.t0
            
            deltaY=obj.eventstress(evtIndex);
            
            % coseismic shear stress change in strike direction
            ynew(3:obj.flt.dgf:(obj.flt.N*obj.flt.dgf))= ...
                ynew(3:obj.flt.dgf:(obj.flt.N*obj.flt.dgf))+deltaY(1:3:(obj.flt.N*3));
            
            % coseismic shear stress change in dip direction
            ynew(4:obj.flt.dgf:(obj.flt.N*obj.flt.dgf))= ...
                ynew(4:obj.flt.dgf:(obj.flt.N*obj.flt.dgf))+deltaY(2:3:obj.flt.N*3);
            
            if isobject(obj.shz)
                % coseismic change in s11
                ynew((obj.flt.N*obj.flt.dgf+ 7):obj.shz.dgf:end)= ...
                    ynew((obj.flt.N*obj.flt.dgf+ 7):obj.shz.dgf:end)+deltaY(1+obj.flt.N*3:6:end);

                % coseismic change in s12
                ynew((obj.flt.N*obj.flt.dgf+ 8):obj.shz.dgf:end)= ...
                    ynew((obj.flt.N*obj.flt.dgf+ 8):obj.shz.dgf:end)+deltaY(2+obj.flt.N*3:6:end);

                % coseismic change in s13
                ynew((obj.flt.N*obj.flt.dgf+ 9):obj.shz.dgf:end)= ...
                    ynew((obj.flt.N*obj.flt.dgf+ 9):obj.shz.dgf:end)+deltaY(3+obj.flt.N*3:6:end);

                % coseismic change in s22
                ynew((obj.flt.N*obj.flt.dgf+10):obj.shz.dgf:end)= ...
                    ynew((obj.flt.N*obj.flt.dgf+10):obj.shz.dgf:end)+deltaY(4+obj.flt.N*3:6:end);

                % coseismic change in s23
                ynew((obj.flt.N*obj.flt.dgf+11):obj.shz.dgf:end)= ...
                    ynew((obj.flt.N*obj.flt.dgf+11):obj.shz.dgf:end)+deltaY(5+obj.flt.N*3:6:end);

                % coseismic change in s33
                ynew((obj.flt.N*obj.flt.dgf+12):obj.shz.dgf:end)= ...
                    ynew((obj.flt.N*obj.flt.dgf+12):obj.shz.dgf:end)+deltaY(6+obj.flt.N*3:6:end);
            end
            % decrease time step after event
            h=min(1e-5,absh);
            absh=abs(h);
            tnew=t+h;
            
            if verbose
                fprintf('coseismic slip distribution ''%s'' at time t=%08.2f yr\n',...
                    obj.evt{evtIndex}.src.name,obj.evt{evtIndex}.src.t0);
            end
            evtIndex=evtIndex+1;
        end
    end
    
    %     if haveEventFcn
    %         [te,ye,ie,valt,stop] = ...
    %             odezero(@ntrp45,eventFcn,eventArgs,valt,t,y,tnew,ynew,t0,h,f,idxNonNegative);
    %         if ~isempty(te)
    %             if output_sol || (nargout > 2)
    %                 teout = [teout, te];
    %                 yeout = [yeout, ye];
    %                 ieout = [ieout, ie];
    %             end
    %             if stop               % Stop on a terminal event.
    %                 % Adjust the interpolation data to [t te(end)].
    %
    %                 % Update the derivatives using the interpolating polynomial.
    %                 taux = t + (te(end) - t)*A;
    %                 [~,f(:,2:7)] = ntrp45(taux,t,y,[],[],h,f,idxNonNegative);
    %
    %                 tnew = te(end);
    %                 ynew = ye(:,end);
    %                 h = tnew - t;
    %                 done = true;
    %             end
    %         end
    %     end
    
    %     if output_sol
    %         nout = nout + 1;
    %         if nout > length(obj.t)
    %             obj.t = [obj.t, zeros(1,chunk,dataType)];  % requires chunk >= refine
    %             obj.y = [obj.y, zeros(neq,chunk,dataType)];
    %             f3d  = cat(3,f3d,zeros(neq,7,chunk,dataType));
    %         end
    %         obj.t(nout) = tnew;
    %         obj.y(:,nout) = ynew;
    %         f3d(:,:,nout) = f;
    %     end
    
    % allocate more memory for max velocity
    if isobject(obj.flt)
        if (1+nMax) > length(obj.flt.tMax)
            obj.flt.tMax=[obj.flt.tMax, zeros(1,chunkPoints,dataType)];
            obj.flt.vMax=[obj.flt.vMax, zeros(1,chunkPoints,dataType)];
        end
    end
       
    if isobject(obj.shz)
        if (1+nMax) > length(obj.shz.tMax)
            obj.shz.tMax=[obj.shz.tMax, zeros(1,chunkPoints,dataType)];
            obj.shz.eMax=[obj.shz.eMax, zeros(1,chunkPoints,dataType)];
        end
    end
    
    % evaluate maximum velocity
    yp=f(:,7);

    if isobject(obj.flt)
        obj.flt.tMax(nMax+1)=tnew;
        vMax=max(sqrt(yp(1:obj.flt.dgf:obj.flt.dgf*obj.flt.N).^2+yp(2:obj.flt.dgf:obj.flt.dgf*obj.flt.N).^2));
        obj.flt.vMax(nMax+1)=vMax;
    end
       
    if isobject(obj.shz)
        obj.shz.tMax(nMax+1)=tnew;
        eMax=max(sqrt(   yp(1+obj.flt.dgf*obj.flt.N:obj.shz.dgf:end).^2 ...
                      +2*yp(2+obj.flt.dgf*obj.flt.N:obj.shz.dgf:end).^2 ...
                      +2*yp(3+obj.flt.dgf*obj.flt.N:obj.shz.dgf:end).^2 ...
                        +yp(4+obj.flt.dgf*obj.flt.N:obj.shz.dgf:end).^2 ...
                      +2*yp(5+obj.flt.dgf*obj.flt.N:obj.shz.dgf:end).^2 ...
                        +yp(6+obj.flt.dgf*obj.flt.N:obj.shz.dgf:end).^2));
        obj.shz.eMax(nMax+1)=eMax;
    end
                       
    if ~isempty(oDir)
       fprintf(timeFID,'%20.16e %20.16e\n',[tnew h]);
       if isobject(obj.flt)
           fprintf(vMaxFID,'%20.16e %20.16e\n',[tnew vMax]);
       end
       if isobject(obj.shz)
           fprintf(eMaxFID,'%20.16e %20.16e\n',[tnew eMax]);
       end
    end

    nMax=nMax+1;
    
    if isobject(obj.flt)
        for k=1:length(obj.flt.observationPoints)
            % allocate more memory for observation points
            if (1+obj.flt.observationPoints{k}.N) > length(obj.flt.observationPoints{k}.t)
                obj.flt.observationPoints{k}.t=[obj.flt.observationPoints{k}.t, zeros(1,chunkPoints,dataType)];
                obj.flt.observationPoints{k}.y=[obj.flt.observationPoints{k}.y, zeros(obj.flt.dgf,chunkPoints,dataType)];
            end
            
            obj.flt.observationPoints{k}.t(obj.flt.observationPoints{k}.N+1)=tnew;
            sample=ynew(1+(obj.flt.observationPoints{k}.index-1)*obj.flt.dgf:obj.flt.observationPoints{k}.index*obj.flt.dgf);
            obj.flt.observationPoints{k}.y(:,obj.flt.observationPoints{k}.N+1)=sample;
            obj.flt.observationPoints{k}.N=obj.flt.observationPoints{k}.N+1;
            
            if ~isempty(oDir)
                fprintf(fltOpts{k},'%20.16e %20.16e',tnew,h);
                fprintf(fltOpts{k},'%21.16e ',sample);
                fprintf(fltOpts{k},'\n');
            end
        end
    end
    
    if isobject(obj.shz)
        for k=1:length(obj.shz.observationPoints)
            % allocate more memory for observation points
            if (1+obj.shz.observationPoints{k}.N) > length(obj.shz.observationPoints{k}.t)
                obj.shz.observationPoints{k}.t=[obj.shz.observationPoints{k}.t, zeros(1,chunkPoints,dataType)];
                obj.shz.observationPoints{k}.y=[obj.shz.observationPoints{k}.y, zeros(obj.shz.dgf,chunkPoints,dataType)];
            end
            
            obj.shz.observationPoints{k}.t(obj.shz.observationPoints{k}.N+1)=tnew;
            sample=ynew(obj.flt.N*obj.flt.dgf+1+(obj.shz.observationPoints{k}.index-1)*obj.shz.dgf: ...
                obj.flt.N*obj.flt.dgf+  (obj.shz.observationPoints{k}.index  )*obj.shz.dgf);
            obj.shz.observationPoints{k}.y(:,obj.shz.observationPoints{k}.N+1)=sample;
            obj.shz.observationPoints{k}.N=obj.shz.observationPoints{k}.N+1;
            
            if ~isempty(oDir)
                fprintf(shzOpts{k},'%20.16e %20.16e',tnew,h);
                fprintf(shzOpts{k},'%22.16e ',sample);
                fprintf(shzOpts{k},'\n');
            end
        end
    end
    
    if isobject(obj.flt)
        for k=1:length(obj.flt.eventCatalogue)
            if (obj.flt.eventCatalogue{k}.vStart<obj.flt.vMax(nMax))
                if (~obj.flt.eventCatalogue{k}.isEvent)
                    % start of a new event
                    obj.flt.eventCatalogue{k}.nEvents=1+obj.flt.eventCatalogue{k}.nEvents;
                    obj.flt.eventCatalogue{k}.evt{obj.flt.eventCatalogue{k}.nEvents}= ...
                        unicycle.ode.evt(tnew,[],nMax,[],ynew(1:obj.flt.dgf*obj.flt.N),[]);
                    obj.flt.eventCatalogue{k}.isEvent=true;
                end
            else if (obj.flt.eventCatalogue{k}.vEnd>obj.flt.vMax(nMax)) && (obj.flt.eventCatalogue{k}.isEvent)
                    % end of current event
                    obj.flt.eventCatalogue{k}.evt{obj.flt.eventCatalogue{k}.nEvents}.tEnd=tnew;
                    obj.flt.eventCatalogue{k}.evt{obj.flt.eventCatalogue{k}.nEvents}.iEnd=nMax;
                    obj.flt.eventCatalogue{k}.evt{obj.flt.eventCatalogue{k}.nEvents}.yEnd= ...
                        ynew(1:obj.flt.dgf*obj.flt.N);
                    obj.flt.eventCatalogue{k}.isEvent=false;
                    fprintf('# catalogue %d, event %d, i=%d, ti=%f, tf=%f\n', ...
                        k, ...
                        obj.flt.eventCatalogue{k}.nEvents, ...
                        obj.flt.eventCatalogue{k}.evt{obj.flt.eventCatalogue{k}.nEvents}.tStart, ...
                        obj.flt.eventCatalogue{k}.evt{obj.flt.eventCatalogue{k}.nEvents}.tEnd);
                end
            end
        end
    end
                    
    if output_ty || haveOutputFcn
        switch outputAt
            case 'Subset'
                nout_new = 0;
                if feval(ioselection,nsteps)
                    nout_new = 1;
                    tout_new = tnew;
                    yout_new = ynew;
                end
            case 'SolverSteps'        % computed points, no refinement
                nout_new = 1;
                tout_new = tnew;
                yout_new = ynew;
                if verbose
                    textprogressbar((tnew-t0)/tspan(end)*100);
                end
            case 'RefinedSteps'       % computed points, with refinement
                tref = t + (tnew-t)*S;
                nout_new = refine;
                tout_new = [tref, tnew];
                yout_new = [ntrp45(tref,t,y,[],[],h,f,idxNonNegative), ynew];
            case 'RequestedPoints'    % output only at tspan points
                nout_new =  0;
                tout_new = [];
                yout_new = [];
                while next <= ntspan
                    if tdir * (tnew - tspan(next)) < 0
                        if haveEventFcn && stop     % output tstop,ystop
                            nout_new = nout_new + 1;
                            tout_new = [tout_new, tnew];
                            yout_new = [yout_new, ynew];
                        end
                        break;
                    end
                    nout_new = nout_new + 1;
                    tout_new = [tout_new, tspan(next)];
                    if tspan(next) == tnew
                        yout_new = [yout_new, ynew];
                    else
                        yout_new = [yout_new, obj.ntrp45(tspan(next),t,y,[],[],h,f,idxNonNegative)];
                    end
                    next = next + 1; 
                end
                if verbose
                    textprogressbar(next/length(tspan)*100);
                end
        end
        
        if nout_new > 0
            if output_ty
                oldnout = nout;
                nout = nout + nout_new;
                if nout > length(obj.t)
                    obj.t = [obj.t, zeros(1,chunk,dataType)];  % requires chunk >= refine
                    obj.y = [obj.y, zeros(neq,chunk,dataType)];
                end
                
                idx = oldnout+1:nout;
                obj.t(idx) = tout_new;
                obj.y(:,idx) = yout_new;
                
            end
            if haveOutputFcn
                stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
                if stop
                    done = true;
                end
            end
        end
    end
    
    if done
        break
    end
    
    % If there were no failures compute a new h.
    if nofailed
        % Note that absh may shrink by 0.8, and that err may be 0.
        temp = 1.25*(err/rtol)^pow;
        if temp > 0.2
            absh = absh / temp;
        else
            absh = 5.0*absh;
        end
    end
    
    % Advance the integration one step.
    t = tnew;
    y = ynew;
    if normcontrol
        normy = normynew;
    end
    if NNreset_f7
        % Used f7 for unperturbed solution to interpolate.
        % Now reset f7 to move along constraint.
        f(:,7) = feval(odeFcn,tnew,ynew,odeArgs{:});
        nfevals = nfevals + 1;
    end
    f(:,1) = f(:,7);  % Already have f(tnew,ynew)
    
end

% free memory
yold=[];

% close output files
if ~isempty(oDir)
    fclose(timeFID);
    if isobject(obj.flt)
        fclose(vMaxFID);
        for k=1:length(obj.flt.observationPoints)
            fclose(fltOpts{k});
        end
        % remove unused memory for observation points
        for k=1:length(obj.flt.observationPoints)
            obj.flt.observationPoints{k}.t=obj.flt.observationPoints{k}.t(1:obj.flt.observationPoints{k}.N);
            obj.flt.observationPoints{k}.y=obj.flt.observationPoints{k}.y(:,1:obj.flt.observationPoints{k}.N);
        end
        % remove unused memory for vMax and tMax
        obj.flt.tMax=obj.flt.tMax(1:nMax);
        obj.flt.vMax=obj.flt.vMax(1:nMax);
    end
    if isobject(obj.shz)
        fclose(eMaxFID);
        for k=1:length(obj.shz.observationPoints)
            fclose(shzOpts{k});
        end
        % remove unused memory for observation points
        for k=1:length(obj.shz.observationPoints)
            obj.shz.observationPoints{k}.t=obj.shz.observationPoints{k}.t(1:obj.shz.observationPoints{k}.N);
            obj.shz.observationPoints{k}.y=obj.shz.observationPoints{k}.y(:,1:obj.shz.observationPoints{k}.N);
        end
        % remove unused memory for eMax and tMax
        obj.shz.tMax=obj.shz.tMax(1:nMax);
        obj.shz.eMax=obj.shz.eMax(1:nMax);
    end
end

solver_output = obj.odefinalize(solver_name, sol,...
    outputFcn, outputArgs,...
    printstats, [nsteps, nfailed, nfevals],...
    nout,...
    haveEventFcn, teout, yeout, ieout,...
    {f3d,idxNonNegative});
if nargout > 0
    varargout = solver_output;
end

if verbose
    textprogressbar('');
    fprintf('\n');
end

end


