function options = odeset(varargin)
%ODESET Create/alter ODE OPTIONS structure.
%   OPTIONS = ODESET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names. 
%   
%   OPTIONS = ODESET(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   OPTIONS = ODESET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties. 
%   
%   ODESET with no input arguments displays all property names and their
%   possible values.
%   
%ODESET PROPERTIES
%   
%RelTol - Relative error tolerance  [ positive scalar {1e-3} ]
%   This scalar applies to all components of the solution vector, and
%   defaults to 1e-3 (0.1% accuracy) in all solvers.  The estimated error in
%   each integration step satisfies e(i) <= max(RelTol*abs(y(i)),AbsTol(i)).
%
%AbsTol - Absolute error tolerance  [ positive scalar or vector {1e-6} ]
%   A scalar tolerance applies to all components of the solution vector.
%   Elements of a vector of tolerances apply to corresponding components of
%   the solution vector. AbsTol defaults to 1e-6 in all solvers. See RelTol.
%   
%NormControl -  Control error relative to norm of solution  [ on | {off} ]
%   Set this property 'on' to request that the solvers control the error in
%   each integration step with norm(e) <= max(RelTol*norm(y),AbsTol). By
%   default the solvers use a more stringent component-wise error control. 
%   
%Refine - Output refinement factor  [ positive integer ]
%   This property increases the number of output points by the specified
%   factor producing smoother output. Refine defaults to 1 in all solvers 
%   except ODE45, where it is 4. Refine does not apply if length(TSPAN) > 2 
%   or the ODE solver returns the solution as a structure.
%   
%OutputFcn - Installable output function  [ function_handle ]
%   This output function is called by the solver after each time step. When
%   a solver is called with no output arguments, OutputFcn defaults to 
%   @odeplot. Otherwise, OutputFcn defaults to [].
%   
%OutputSel - Output selection indices  [ vector of integers ]
%   This vector of indices specifies which components of the solution vector
%   are passed to the OutputFcn. OutputSel defaults to all components.
%   
%Stats - Display computational cost statistics  [ on | {off} ]
%   
%Jacobian - Jacobian function [ function_handle | constant matrix ]
%   Set this property to @FJac if FJac(t,y) returns dF/dy, or to
%   the constant value of dF/dy.   
%   For ODE15I solving F(t,y,y') = 0, set this property to @FJac if
%   [dFdy, dFdyp] = FJac(t,y,yp), or to a cell array of constant
%   values {dF/dy,dF/dyp}.
%      
%JPattern - Jacobian sparsity pattern [ sparse matrix ]
%   Set this property to a sparse matrix S with S(i,j) = 1 if component i of
%   F(t,y) depends on component j of y, and 0 otherwise.
%   For ODE15I solving F(t,y,y') = 0, set this property to
%   {dFdyPattern,dFdypPattern}, the sparsity patterns of dF/dy and
%   dF/dy', respectively. 
%   
%Vectorized - Vectorized ODE function  [ on | {off} ]
%   Set this property 'on' if the ODE function F is coded so that 
%   F(t,[y1 y2 ...]) returns [F(t,y1) F(t,y2) ...]. 
%   For ODE15I solving F(t,y,y') = 0, set this property to
%   {yVect,ypVect}. Setting yVect 'on' indicates that 
%   F(t,[y1 y2 ...],yp) returns [F(t,y1,yp) F(t,y2,yp) ...].  
%   Setting ypVect 'on' indicates that F(t,y,[yp1 yp2 ...])
%   returns [F(t,y,yp1) F(t,y,yp2) ...].   
%      
%Events - Locate events  [ function_handle ]
%   To detect events, set this property to the event function.
%   
%Mass - Mass matrix [ constant matrix | function_handle ]
%   For problems M*y' = f(t,y) set this property to the value of the constant
%   mass matrix. For problems with time- or state-dependent mass matrices,
%   set this property to a function that evaluates the mass matrix.
%
%MStateDependence - Dependence of the mass matrix on y [ none | {weak} | strong ] 
%   Set this property to 'none' for problems M(t)*y' = F(t,y). Both 'weak' and
%   'strong' indicate M(t,y), but 'weak' will result in implicit solvers
%   using approximations when solving algebraic equations.
%   
%MassSingular - Mass matrix is singular  [ yes | no | {maybe} ]
%   Set this property to 'no' if the mass matrix is not singular.
%
%MvPattern - dMv/dy sparsity pattern [ sparse matrix ]
%   Set this property to a sparse matrix S with S(i,j) = 1 if for any k, the 
%   (i,k) component of M(t,y) depends on component j of y, and 0 otherwise.  
%
%InitialSlope - Consistent initial slope yp0 [ vector ]
%   yp0 satisfies M(t0,y0)*yp0 = F(t0,y0).
%
%InitialStep - Suggested initial step size  [ positive scalar ]
%   The solver will try this first.  By default the solvers determine an
%   initial step size automatically. 
%   
%MaxStep - Upper bound on step size  [ positive scalar ]
%   MaxStep defaults to one-tenth of the tspan interval in all solvers.
%
%NonNegative - Non-negative solution components [ vector of integers ]
%   This vector of indices specifies which components of the
%   solution vector must be non-negative.  NonNegative defaults to []. 
%   This property is not available in ODE23S, ODE15I.  In ODE15S,
%   ODE23T, and ODE23TB, the property is not available for problems
%   where there is a mass matrix.
%
%BDF - Use Backward Differentiation Formulas in ODE15S  [ on | {off} ]
%   This property specifies whether the Backward Differentiation Formulas
%   (Gear's methods) are to be used in ODE15S instead of the default
%   Numerical Differentiation Formulas. 
%
%Verbose - display a progress bar during simulation [ {true} | false ]
%
%MaxOrder - Maximum order of ODE15S and ODE15I [ 1 | 2 | 3 | 4 | {5} ]
%   
%   See also ODEGET, ODE45, ODE23, ODE113, ODE15I, ODE15S, ODE23S, ODE23T, ODE23TB,
%            FUNCTION_HANDLE.

%   Mark W. Reichelt and Lawrence F. Shampine, 5/6/94
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.46.4.15 $  $Date: 2011/05/17 02:23:18 $

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
  fprintf('          AbsTol: [ positive scalar or vector {1e-6} ]\n');
  fprintf('          RelTol: [ positive scalar {1e-3} ]\n');
  fprintf('     NormControl: [ on | {off} ]\n');
  fprintf('     NonNegative: [ vector of integers ]\n'); 
  fprintf('            oDir: [ output directory ]\n'); 
  fprintf('       OutputFcn: [ function_handle ]\n');
  fprintf('       OutputSel: [ vector of integers ]\n');
  fprintf('          Refine: [ positive integer ]\n');  
  fprintf('           Stats: [ on | {off} ]\n');
  fprintf('     InitialStep: [ positive scalar ]\n');
  fprintf('         MaxStep: [ positive scalar ]\n');
  fprintf('             BDF: [ on | {off} ]\n');
  fprintf('         Verbose: [ {true} | false ]\n');
  fprintf('        MaxOrder: [ 1 | 2 | 3 | 4 | {5} ]\n');
  fprintf('        Jacobian: [ matrix | function_handle ]\n');
  fprintf('        JPattern: [ sparse matrix ]\n');
  fprintf('      Vectorized: [ on | {off} ]\n');
  fprintf('            Mass: [ matrix | function_handle ]\n');
  fprintf('MStateDependence: [ none | {weak} | strong ]\n');
  fprintf('       MvPattern: [ sparse matrix ]\n');
  fprintf('    MassSingular: [ yes | no | {maybe} ]\n');
  fprintf('    InitialSlope: [ vector ]\n');
  fprintf('          Events: [ function_handle ]\n');
  fprintf('          Netout: [NetCDF output filename ] \n');
  fprintf('\n');
  return;
end

Names = [
    'AbsTol          '
    'BDF             '
    'Verbose         '
    'Events          '
    'InitialStep     '
    'Jacobian        '
    'JConstant       '             % backward compatibility
    'JPattern        '
    'Mass            '
    'MassSingular    '
    'MaxOrder        '
    'MaxStep         '
    'NonNegative     ' 
    'NormControl     '
    'oDir            '
    'OutputFcn       '
    'OutputSel       '
    'Refine          '
    'RelTol          '
    'Stats           '
    'Vectorized      '
    'MStateDependence'
    'MvPattern       '
    'InitialSlope    '
    'Netout          '
    ];
m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(message('MATLAB:odeset:NoPropNameOrStruct', i));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error(message('MATLAB:odeset:ArgNameValueMismatch'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error(message('MATLAB:odeset:NoPropName', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(message('MATLAB:odeset:InvalidPropName', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
            matches = deblank(Names(j(1),:));
        for k = j(2:length(j))'
                matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
        end
            error(message('MATLAB:odeset:AmbiguousPropName',arg,matches));
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(message('MATLAB:odeset:NoValueForProp', arg));
end
