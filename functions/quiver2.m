function [HQ,HCB,HL] = quiver2(varargin)
%QUIVER   Same as quiver but with customizations.
%
%   SYNTAX:
%                   quiver2(U,V)
%                   quiver2(U,V,W)
%                   quiver2(U,V,S)
%                   quiver2(U,V,W,S)
%                   quiver2(x,y,U,V)
%                   quiver2(x,y,z,U,V,W)
%                   quiver2(x,y,z,U,V,W,S)
%                   quiver2(...,LINESPEC)
%                   quiver2(...,'filled')
%                   quiver2(...,'c=',CMAP)        % colored arrows
%                   quiver2(...,'w=',WLIM)        % widened
%                   quiver2(...,'n=',MAG)         % normalized
%                   quiver2(...,'l=',LLIM)        % limitized
%                   quiver2(...,'s=','screen')    % slanted by (U,V) angle!
%                   quiver2(...,'t=','quiver')    % try 'feather' or 'b'
%                   quiver2(...,'a@','fancy')     % Uses ARROW function!
%                   quiver2(...,'x@',SAXIS)       % try 'x'
%                   quiver2(...,'f@','Edge')      % try 'Face'
%                   quiver2(...,'l@',Length)      % see ARROW for details
%                   quiver2(...,'b@',BaseAngle)
%                   quiver2(...,'t@',TipAngle)
%                   quiver2(...,'w@',Width)
%                   quiver2(...,'p@',Page)
%                   quiver2(...,'c@',CrossDir)
%                   quiver2(...,'n@',NormalDir)
%                   quiver2(...,'e@',Ends)
%                   quiver2(...,P/V) % QUIVER optional pairs input
%                   quiver2(AX,...)
%              HQ = quiver2(...);
%        [HQ,HCB] = quiver2(...);
%     [HQ,HCB,HL] = quiver2(...);
%
%   INPUT:
%     U        - X component (with/without NaNs).
%     V        - Y component (with/without NaNs).
%     W        - Z component (with/without NaNs).
%                DEFAULT: [] (not used)
%     x        - X axis (matrix or vector).
%                DEFAULT: 1:size(U,2)
%     y        - Y axis (matrix or vector).
%                DEFAULT: 1:size(U,1)
%     z        - Z axis (matrix or vector).
%                DEFAULT: 1
%     S        - Arrow scaling (see QUIVER for details). Recommended not to
%                use it and use 'n=',MAG instead.
%                DEFAULT: 0 (do not scales)
%     LINESPEC - Arrow line and marker style (see QUIVER and PLOT for
%                details). For example: '.' draws lines instead of arrows.
%                DEFAULT: (draws normal arrows)
%     'filled' - Fills the arrow markers.
%                DEFAULT (do not fills)
%     'c='     - Colors arrows according to (U,V) magnitude with the
%                colormap CMAP. CMAP may be [] to use blue color. 
%                DEFAULT: (is used by default with the current colormap)
%     'w='     - Widens arrows according to (U,V) magnitude with
%                WLIM=[Wmin Wmax] pixel limits widths.  
%                DEFAULT: (is used by default with WLIM=[0.5 2.5])
%     'n='     - Draws arrows with normalized lengths to MAG. Use MAG=0 to
%                not normalize.
%                DEFAULT: 1 (by default the arrows are normalized!)
%     'l='     - Draws arrows squeezed to this LLIM limits for length.
%                DEFAULT: (not used by default)
%     's='     - Slants arrows accordingly to (U,V) angle. That is, a 45°
%                arrow will do have a 45° slope on 'screen'. Use 'normal'
%                for default QUIVER slant.
%                DEFAULT: (is used by default with 'screen')
%     't='     - Function type. One of 'quiver' or 'feather' way (adds an
%                horizontal line). May be a LINESPEC instead of 'feather'!.
%                Use 'feather-date' or [LINESPEC '-date'] to treat xticks
%                as dates.
%                DEFAULT: (used with 'quiver' but 'feather' if x or y = [])
%     'a@'     - ARROWs custom type. One of 'fancy', '90' or 'normal'.
%                Write one of yours inside this file. See ARROW for details
%                in this '@' optional inputs.
%                DEFAULT: (not used)
%     'x@'     - ARROWs with tips parallel to specified axes: 'x', 'y' or
%                'z'. Or negative as '-x'.
%                DEFAULT: (not used)
%     'f@'     - ARROWs with colored 'face' (default) or 'edge'. Use of
%                'filled' colored both.
%                DEFAULT: (not used)
%     'l@'     - ARROWs Length property.
%                DEFAULT: (not used)
%     'b@'     - ARROWs BaseLine property.
%                DEFAULT: (not used)
%     't@'     - ARROWs TipAngle property.
%                DEFAULT: (not used)
%     'w@'     - ARROWs Width property.
%                DEFAULT: (not used)
%     'p@'     - ARROWs Page property. 
%                DEFAULT: (not used)
%     'c@'     - ARROWs CrossDir property. 
%                DEFAULT: (not used)
%     'n@'     - ARROWs NormalDir property.
%                DEFAULT: (not used)
%     'e@'     - ARROWs Ends property.
%                DEFAULT: (not used)
%     P/V      - Property/Value optional QUIVER function pairs inputs.
%                DEFAULT: (not use any)
%     AX       - Draws arrows on axes with handle AX instead of current
%                one. 
%                DEFAULT: gca (current axes)
%
%   OUTPUT:
%     HQ  - Returns quiver handle.
%     HCB - Returns generated COLORBAR handle.
%     HL  - Handle of line when 'feather' is used.
%
%   DESCRIPTION:
%     This programs is the same as QUIVER and QUIVER3 but with the colors,
%     lengths, widths and arrows customization.
%
%     Besides, it works like FEATHER function when y is [].
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * ADDITIONAL NOTES are included inside this file.
%
%   EXAMPLE:
%     % EXAMPLE 1. quiver2 vs quiver
%      % Data (45° field + NaNs):
%       Nx=50; Ny=10; PNAN=50;
%       [x,y]   = meshgrid((1:Nx)/2,(1:Ny)*2);
%       U       = 10*reshape(sort(rand(Nx*Ny,1)),Ny,Nx);
%       A       = repmat(45*pi/180,Ny,Nx);
%       [u,v]   = pol2cart(A,U); 
%       inan    = randperm(Nx*Ny); inan = inan(1:floor((PNAN/100)*Nx*Ny));
%       u(inan)=NaN; v(inan)=NaN;
%      % quiver2:
%       a = subplot(221);
%       quiver2(x,y,u,v,'filled','n=',2)   % augmented length by 2
%       title '45° colored QUIVER2 field vs ...'
%      % quiver:
%       b = subplot(223);
%       quiver(x,y,u,v), title '... same? 45° QUIVER field'
%       linkaxes([a b]), zoom on
%
%     % EXAMPLE 2. quiver2 vs feather
%       theta=(-90:10:90)*pi/180; r=2*ones(size(theta)).*theta;
%       t     = (0:(length(theta)-1))/24 + now; % time/date axis
%       [u,v] = pol2cart(theta,r);
%       subplot(222)
%        set(gca,'XTick',t(1:4:end))
%        quiver2(t,[],u,v,'t=','k-date','x@','x')
%        title 'QUIVER2 + ARROW vs ...'
%       subplot(224)
%        feather(u,v), title '... FEATHER'
%     
%   SEE ALSO:
%     QUIVER, QUIVER3, FEATHER
%     and
%     ARROW by Erik Johnson and CMAPPING, TLABEL by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   quiver2.m
%   VERSION: 1.2 (Nov 12, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   ADDITIONAL NOTES:
%     * By default, QUIVER2 normalizes, colored and widens arrows
%       according to its length. To get normal behavior of QUIVER or
%       QUIVER3 use empties on '=' entries or 
%         quiver2(...,'c=','b','w=',0.5,'s=','normal')
%       and 't=','feather' respectively.
%     * Default properties of QUIVER2 are at the top of this file, type
%         edit quiver2
%       to see/change them. For example, if ARROW is used and its 'Length'
%       property not ('l@'), this programs generates tip arrows with 33%
%       the size of the whole arrow (like QUIVER does). There, you can
%       change it to 50% with alpha=0.50.
%     * ARROW function by Erik Johnson is highly recommended on the
%       Mathworks FileExchange for arrows creation. Download it and make
%       your own customized arrows with the '@' optional inputs. Type 
%         arrow properties
%       to get more details.
%     * This programs may also works with CMAPPING and TLABEL from the
%       FileExchange.

%   REVISIONS:
%   1.0      Released. (Jun 30, 2009)
%   1.1      Fixed bug with inputs check in and on a error message, thanks
%            to Ayal Anis. (Oct 23, 2009)
%   1.2      Fixed bug with MAG=0 value, thanks to Ayal Anis. Spell
%            checked. (Nov 12, 2009) 

%   DISCLAIMER:
%   quiver2.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Defaults.
AX        = [];        % Gets current axes.
CMAP      = [];        % Gets figure colormap.
WLIM      = [0.5 2.5]; % Line width limits according to U,V length.
MAG       = 1;         % Normalization.
LLIM      = [];        % Not used by default.
SLANT     = 'screen';  % Slant 'screen' or 'normal'
TYPE      = 'quiver';  % May change according to x,y inputs.
S         = 0;         % Scaling factor
LINESPEC  = {};        % LINESPEC for line arrows
filled    = {};        % 'filled'?
visible   = 'on';      % 'visible'?
K         = 0.1;       % Increase factor for limits
Aopt      = {};        % ARROWs default P/V options ------
alpha     = 0.33;      % Length of arrow head with respect of length
EorF      = 'Edge';    % Draw 'Face' or 'Edge'
CUSTOM    = [];        % 'fancy', for example
SAXES     = [];        % 'x' for example

% Customized ARROWs (add your owns).
Anames = {'fancy' ,{'BaseAngle',50};...
          '90'    ,{'BaseAngle',90};...
          'normal',{'Length'   ,[],...
                    'BaseAngle',[],...
                    'TipAngle' ,[],...
                    'Width'    ,[],...
                    'Page'     ,[],...
                    'CrossDir' ,[],...
                    'NormalDir',[],...
                    'Ends'     ,[]}};

% Checks number of inputs and outputs.
if     nargin<2
 error('CVARGAS:quiver2:notEnoughInputs',...
  'At least 2 inputs are required.')
elseif nargout>3 
 error('CVARGAS:quiver2:tooManyOutputs',...
  'At most 3 output are allowed.')
end

% Checks AX input.
if (length(varargin{1})==1) && ~ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'Type'),'axes')
 AX = varargin{1}; varargin(1) = [];
 if length(varargin)<2
  error('CVARGAS:quiver2:notEnoughInputsAfterHandle',...
  'At least 2 inputs plus handle are required.')
 end
end

% Checks numeric x,y,z,U,V,W,S 7 inputs.
% Note: this is done here to reduce memory use.
% Fixed bug, thanks to Ayal Anis. Oct 2009.
Nopt = {[],[],[],[],[],[],S};
temp = 0;
for k = 1:min(7,length(varargin)) 
 if ischar(varargin{k}), temp = 1; break, end 
end
k = k-temp;
if k<2
 error('CVARGAS:quiver2:incorrectNumericInputs',...
  'At least the first two inputs must be numerical.')
end
if k==2,                    Nopt(4:5) = varargin(1:2);
elseif k==3,                Nopt(4:5) = varargin(1:2);
 if numel(varargin{3})~=1,  Nopt(6)   = varargin(3);
 else                       Nopt(7)   = varargin(3);
  if numel(varargin{1})==1
   error('CVARGAS:quiver2:undefinedNumericOption',...
    ['3rd numerical option may be W or S. '...
    'Use empties x,y (and z, W) to avoid confusion.'])
  end
 end
elseif k==4
 if numel(varargin{4})==1,  Nopt(4:7) = varargin(1:4);
  if numel(varargin{3})==1
   error('CVARGAS:quiver2:undefinedNumericOption',...
    ['4rd numerical option may be V or S. '...
    'Use empties x,y (and z, W) to avoid confusion.'])
  end
 else                       Nopt([1 2 4 5])   = varargin(1:4);
 end
elseif k==5,                Nopt([1 2 4 5 7]) = varargin(1:5);
elseif k==6,                Nopt(1:6)         = varargin(1:6);
elseif k==7,                Nopt(1:7)         = varargin(1:7);
end
varargin(1:k) = [];
[x,y,z,U,V,W,S] = Nopt{:}; clear Nopt

% Checks U,V,W.
su=size(U);  sv=size(V);  sw=size(W);
nu=prod(su); nv=prod(sv); nw=prod(sw);
sm=max([su;sv;sw],[],1);  nm=prod(sm);
if any([length(su) length(sv) length(sw)]~=2)
 error('CVARGAS:quiver2:incorrectUvwDimensions',...
  'U,V(,W) must be 2-dimensional matrix or vectors.')
elseif isempty(U)
 if nargout>0
  HQ=[]; HCB=[]; HL=[];
 end
 return
elseif (nm>1) && any([nu nv nw]==1)
 if nu==1
  U = repmat(U,sm); su=sm; nu=nm;
 elseif nv==1
  V = repmat(V,sm);
 elseif nw==1
  W = repmat(W,sm);
 end
elseif (nu~=nv) || (~isempty(W) && (nu~=nw))
 error('CVARGAS:quiver2:incorrectUvwLengths',...
  'U,V(,W) must be have the same number of elements or W=[].')
else
 if ~all(su==sv)
  su = [nu 1]; % If vector, forces column wisee
  U=U(:); V=V(:);
  if ~isempty(W)
   W=W(:);
  end
 else
  if ~isempty(W) && ~all(su==sw)
   error('CVARGAS:quiver2:incorrectUvwShapes',...
  'U,V,W must be have equal shape.')
  end
 end
end

% Checks S.
if (S<0) || ~isfinite(S)
 error('CVARGAS:quiver2:incorrectSValue',...
  'S must be a positive scalar.')
end
 
% Checks string inputs.
N = length(varargin);
k = 1;
while k<=N
 if isempty(varargin{k}) || ~ischar(varargin{k})
  % Empties!
  error('CVARGAS:quiver2:invalidNotCharInput',...
   'Input (name) must be an string.')
 end
 if (k~=N) && (length(varargin{k})==2) && strcmp(varargin{k}(2),'=')
  % QUIVER2 options with '='.
  switch lower(varargin{k}(1))
   case 'c', CMAP = varargin{k+1}; if isempty(CMAP), CMAP = [0 0 1]; end
   case 'w', WLIM = varargin{k+1}; if isempty(WLIM), WLIM = 0.5; end
   case 'n', MAG  = varargin{k+1}; if isempty(MAG),  MAG  = 0; end
   case 'l', LLIM = varargin{k+1}; if isempty(LLIM), LLIM = []; end
   case 't', TYPE = varargin{k+1}; if isempty(TYPE), TYPE = 'quiver'; end
   case 's', SLANT = varargin{k+1};if isempty(SLANT), SLANT = 'normal';end
   otherwise
    error('CVARGAS:quiver2:incorrectOptionalInputLength',...
    'Optional ''='' input must be preceded by one of ''cwnlts''.')
  end
  varargin(k:k+1) = []; N = N-2; k = k-1;
 elseif (k~=N) && (length(varargin{k})==2) && strcmp(varargin{k}(2),'@')
  % ARROW option with '@'.
  switch lower(varargin{k}(1))
   case 'a', CUSTOM = varargin{k+1};
   case 'x', SAXES = varargin{k+1};
   case 'f', EorF = varargin{k+1};
   case 'l', Aopt = {Aopt{:},'Length',varargin{k+1}};
   case 'b', Aopt = {Aopt{:},'BaseAngle',varargin{k+1}};
   case 't', Aopt = {Aopt{:},'TipAngle',varargin{k+1}};
   case 'w', Aopt = {Aopt{:},'Width',varargin{k+1}};
   case 'p', Aopt = {Aopt{:},'Page',varargin{k+1}};
   case 'c', Aopt = {Aopt{:},'CrossDir',varargin{k+1}};
   case 'n', Aopt = {Aopt{:},'NormalDir',varargin{k+1}};
   case 'e', Aopt = {Aopt{:},'Ends',varargin{k+1}};
   otherwise
    error('CVARGAS:quiver2:incorrectOptionalArrowInputLength',...
    'Optional ''@'' input must be preceded by one of ''axflbtwpcne''.')
  end
  varargin(k:k+1)=[]; N=N-2; k=k-1;
 elseif (k<=2) && strncmpi('filled',varargin{k},length(varargin{k}))
  % QUIVER 'filled' option.
  filled = {'filled'}; 
  varargin(k)=[]; N=N-1; k=k-1;
 elseif (k~=N) && ...
                 strncmpi('visible',varargin{k},max(2,length(varargin{k})))
  % QUIVER P/V option.
  visible = varargin{k+1};
  varargin(k:k+1)=[]; N=N-2; k=k-1;
 elseif (k~=N) && strncmpi('Parent',varargin{k},length(varargin{k}))
  % QUIVER P/V option.
  AX = varargin{k+1};
  varargin(k:k+1)=[]; N=N-2; k=k-1;
 elseif k==1
  % QUIVER LINESPEC option?
  tempf = get(0,'CurrentFigure');
  if ~isempty(tempf), tempa = get(tempf,'CurrentAxes'); end
  hf2 = figure('Visible','off','HandleVisibility','off');
  ha2 = axes('HandleVisibility','off','Parent',hf2);
  try
   hl2 = plot([0 1],[0 1],varargin{k},'Parent',ha2);
   % If no error continues to get LineStyle, Marker and Color.
   LINESPEC = {varargin{k}};
   LineStyle = get(hl2,'LineStyle');
   if strcmp(LineStyle,'-') && isempty(regexp(LINESPEC,'[-]','once')) 
    LineStyle = []; 
   end
   if strcmp(LineStyle,'none'), LineStyle = []; end
   Marker = get(hl2,'Marker');
   if strcmp(Marker,'none'), Marker = []; end
   icol  = regexp(lower(LINESPEC{1}),'[ymcrgbwk]');
   if ~isempty(icol), Color = 'ymcrgbwk'; Color = Color(icol);
   else Color = 'b'; end
   varargin(k)=[]; N=N-1; k=k-1;
  end
  delete(hf2)
  set(0,'CurrentFigure',tempf)
  if ~isempty(tempf), set(tempf,'CurrentAxes',tempa), end
 else
  error('CVARGAS:quiver2:incorrectStringInput',...
   'Incorrect string option.')
 end
 k = k+1;
end

% Checks AX.
if isempty(AX), AX = gca; end

% Checks CMAP.
if isempty(CMAP)
 % LINESPEC given, looks for colors.
 if ~isempty(LINESPEC) && ~isempty(Color)
  switch Color
   case 'y', CMAP = [1 1 0];
   case 'm', CMAP = [1 0 1];
   case 'c', CMAP = [0 1 1];
   case 'r', CMAP = [1 0 0];
   case 'g', CMAP = [0 1 0];
   case 'b', CMAP = [0 0 1];
   case 'w', CMAP = [1 1 1];
   case 'k', CMAP = [0 0 0];
  end
 else CMAP = colormap(AX); 
 end
elseif ~isnumeric(CMAP)
 if exist('cmapping.m','file')==2, CMAP = cmapping([],CMAP);
 else
  warning('CVARGAS:quiver2:cmappingNotFound',...
   {['Not found CMAPPING funtion, get it at ' ...
    'http://www.mathworks.com/matlabcentral/fileexchange']; ...
    'Used current COLORMAP instead of specified one.'})
  CMAP = colormap(AX);
 end
elseif (ndims(CMAP)~=2) || (size(CMAP,2)~=3)
 error('CVARGAS:quiver2:incorrectCmapRgb',...
  'Numeric CMAP must be a valid RGB colormap with 3 columns.')
end

% Checks WLIM.
if length(WLIM)==1, WLIM = [WLIM WLIM]; end
if (length(WLIM)~=2) || ~all(isfinite(WLIM)) || (WLIM(1)>WLIM(2)) || ...
  (WLIM(1)<=0) || (WLIM(2)<=0)
 error('CVARGAS:quiver2:incorrectWlimLimits',...
  'WLIM must have 2 increasing or a single positive scalar(s).')
end

% Checks MAG.
if (length(MAG)~=1) || ~isfinite(MAG) || (MAG<0) % Fixed bug v1.2
 error('CVARGAS:quiver2:incorrectMagValue',...
  'MAG must be a finite not negative scalar.')
end

% Checks LLIM.
if isempty(LLIM), % continue
elseif (length(LLIM)~=2) || (LLIM(1)>LLIM(2)) || (LLIM(1)<0) || ...
  (LLIM(2)==0)
 error('CVARGAS:quiver2:incorrectUvlimValue',...
  'LLIM must be empty [], or 2 positive valid length limits.')
end

% Checks TYPE.
if strncmpi('quiver',TYPE,length(TYPE)), TYPE = 'quiver';
else
 % continue, hopping it is a valid LINESPEC or 'feather' (plus '-d')
end

% Checks SLANT.
if     strncmpi('screen',SLANT,length(SLANT)), SLANT = 'screen';
elseif strncmpi('normal',SLANT,length(SLANT)), SLANT = 'normal';
else
 error('CVARGAS:quiver2:incorrectSlantValue',...
  'SLANT value must be one of ''screen'' or ''normal''.')
end

% Checks x,y,z sizes.
sx=size(x);  sy=size(y);  sz=size(z);
nx=prod(sx); ny=prod(sy); nz=prod(sz);
if any([length(sx) length(sy) length(sz)]~=2)
 % x,y,z must be vectors or 2D-matrix.
 error('CVARGAS:quiver2:incorrectXyDimensions',...
  'x,y(,z) must be 2-dimensional matrix or vectors or [].')
elseif ~all([isempty(W) isempty(z)]==true)
 % W and z must be both empties, not just one.
 error('CVARGAS:quiver2:incorrectWzempty',...
  'Both, W and z must be empty, not just one.')
elseif (nx==nu) && (ny==nu) && (isempty(z) || (nz==nu))
 % Same elements as U,V.
 x = reshape(x,su); y = reshape(y,su);
 if ~isempty(W), z = reshape(z,su); end
elseif (nx==su(2)) && (ny==su(1)) && (isempty(z) || (nz==1))
 % Axis to matrix.
 [x,y] = meshgrid(x(:),y(:));
 if ~isempty(z), z = repmat(z,su); end
elseif isempty(x) && isempty(y) && (isempty(z) || (nz==1))
 % x,y both [].
 [x,y] = meshgrid(1:su(2),1:su(1)); % default grid
 if ~isempty(z), z = repmat(z,su); end
elseif isempty(x) && (ny==nu) && (isempty(z) || (nz==1))
 % x [], but y matrix.
 x = repmat(1:su(2),su(1),1); y = reshape(y,su);
 if ~isempty(z), z = repmat(z,su); end
 if strcmp('quiver',TYPE), TYPE = 'feather'; end
elseif isempty(x) && (ny==su(1)) && (isempty(z) || (nz==1))
 % x [], but y axis.
 [x,y] = meshgrid(ones(su(2),1),y(:));
 if ~isempty(z), z = repmat(z,su); end
 if strcmp('quiver',TYPE), TYPE = 'feather'; end
elseif (nx==nu) && isempty(y) && (isempty(z) || (nz==1))
 % x matrix, but y [].
 x = reshape(x,su); y = repmat((0:su(1)-1)',1,su(2));
 if ~isempty(z), z = repmat(z,su); end
 if strcmp('quiver',TYPE), TYPE = 'feather'; end
elseif (nx==su(2)) && isempty(y) && (isempty(z) || (nz==1))
 % x vector, but y [].
 [x,y] = meshgrid(x(:),zeros(su(1),1));
 if ~isempty(z), z = repmat(z,su); end
 if strcmp('quiver',TYPE), TYPE = 'feather'; end
else
 error('CVARGAS:quiver2:incorrectXyInput',...
  'Incorrect x,y(,z) inputs shapes.')
end

% Checks CUSTOM arrows.
if ~isempty(CUSTOM)
 k = 1;
 while k<=size(Anames,1) 
  if strncmpi(Anames{k,1},CUSTOM,min(length(Anames{k,1}),length(CUSTOM)))
   temp = Anames{k,2}; Aopt = {Aopt{:} temp{:}}; break
  end
  k = k+1;
 end
 if k>size(Anames,1)
  error('CVARGAS:quiver2:unrecognizedCustomArrowsName',...
   ['Unrecognized custom ARROWs name. Must be one of '...
   '''fancy'',''90'',''normal'' or user created inside this function.'])
 end
end
if ~isempty(SAXES)
 switch lower(SAXES)
  case 'x',  vec = [1 0 0];
  case 'y',  vec = [0 1 0];
  case 'z',  vec = [0 0 1];
  case '-x', vec = -[1 0 0];
  case '-y', vec = -[0 1 0];
  case '-z', vec = -[0 0 1];
  otherwise
   error('CVARGAS:quiver2:unrecognizedSaxesArrowsInput',...
   ['Unrecognized SAXES string. Must be one of ''xyz'' or '...
    '''-x'', ''-y'', ''-z''.'])
 end
 Aopt = {Aopt{:},'CrossDir',vec};
end
if strncmpi('face',EorF,length(EorF)),     EorF = 'Face';
elseif strncmpi('edge',EorF,length(EorF)), EorF = 'Edge';
else
 error('CVARGAS:quiver2:unrecognizedFaceEdgeInput',...
  '@f input mus be followed by ''face'' or ''edge''.')
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Checks 2D quiver.
do2d = isempty(W);

% Checks ARROWs.
doarrows = ~isempty(Aopt);
if doarrows && ~(exist('arrow.m','file')==2)
 warning('CVARGAS:quiver2:notFoundArrowFunction',...
  ['You must download ARROW function from the FileExchange '... % Fixed bug
   'to use the ''@'' properties. Normal QUIVER was used instead.'])
  doarrows = false;
end

% Checks complex.
doreal = false;
if ~isreal(U), U = real(U); doreal = true; end
if ~isreal(V), V = real(V); doreal = true; end
if ~isreal(W), W = real(W); doreal = true; end
if ~isreal(x), x = real(x); doreal = true; end
if ~isreal(y), y = real(y); doreal = true; end
if ~isreal(z), z = real(z); doreal = true; end
if doreal
 error('CVARGAS:quiver2:ignoredImaginaryPart',...
  'Ignored imaginary part(s).')
end

% Arrows length.
L = hypot(U,V); if ~do2d, L = hypot(L,W); end

% Do not plots NaN (not even a single dot).
inan = ~isfinite(U) | ~isfinite(V);
if ~do2d, inan = inan | ~isfinite(W); end
x(inan)=NaN; y(inan)=NaN; z(inan)=NaN;
U(inan)=NaN; V(inan)=NaN; W(inan)=NaN; L(inan)=NaN;
if all(inan(:))
 if nargout>0, HQ=[]; HCB=[]; HL=[]; end
 return
end

% Normalizes.
if (MAG~=0) && strcmp('quiver',TYPE)
 bad = (L==0) | inan;
 U(~bad)=U(~bad)./L(~bad)*MAG;     V(~bad)=V(~bad)./L(~bad)*MAG;
 if ~do2d, W(~bad)=W(~bad)./L(~bad)*MAG; end
end

% Gets limits.
Lmin = min(L(:)); Lmax = max(L(:));
if ~isempty(LLIM)
 if isfinite(LLIM(1)), Lmin = LLIM(1); end
 if isfinite(LLIM(2)), Lmax = LLIM(2); end
end

% Color bands.
Nc = size(CMAP,1);
dL = (Lmax-Lmin)/Nc;
if dL<eps*10^6; dL = 0; Nc = 1; end

% Gets modes before any changes.
autoxlim = strcmp(get(AX,'XLimMode' ),'auto');
autoylim = strcmp(get(AX,'YLimMode' ),'auto');
autozlim = strcmp(get(AX,'ZLimMode' ),'auto');
autoxtic = strcmp(get(AX,'XTickMode'),'auto');
autoytic = strcmp(get(AX,'YTickMode'),'auto');
autoztic = strcmp(get(AX,'ZTickMode'),'auto');
autoclim = strcmp(get(AX,'CLimMode' ),'auto');
ihold    = ishold(AX);
if ~autoxlim, limx2  = get(AX,'XLim' ); limx = limx2; end
if ~autoylim, limy2  = get(AX,'YLim' ); limy = limy2; end
if ~autozlim, limz2  = get(AX,'ZLim' ); limz = limz2; end
if ~autoxtic, xtick2 = get(AX,'XTick'); end

% Prepares for x-axis for 'feather'.
if ~strcmp(TYPE,'quiver')
 xlim = [min(x(:)) max(x(:))];
 if ((su(2)==1) && (sum(~inan(:))==1) && all(diff(x(~inan(:)))==0)) || ...
   (xlim(1)==xlim(2))
  x(~inan(:)) = 1;
  if ~autoxlim, limx  = limx2-xlim(1)+1; end
  if ~autoxtic, xtick = xtick2-xlim(1)+1; end
 else
  x(~inan(:)) = interp1(xlim,[1 su(2)],x(~inan(:)),'linear','extrap');
  if ~autoxlim, limx = interp1(xlim,[1 su(2)],limx2 ,'linear','extrap');end
  if ~autoxtic, xtick= interp1(xlim,[1 su(2)],xtick2,'linear','extrap');end
 end
end

% Start to draw arrows. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializes handles.
HQ = NaN(Nc,1);

% First color band.
if Nc==1, ind = (L<Inf);
else      ind = (L<(Lmin+dL)); end
if ~isempty(ind)
 if do2d
  HQ(1) = quiver(AX,x(ind),y(ind),U(ind),V(ind),S,...
      LINESPEC{:},filled{:},varargin{:},'Color',CMAP(1,:),'Visible','off');
 else
  HQ(1) =quiver3(AX,x(ind),y(ind),z(ind),U(ind),V(ind),W(ind),S,...
      LINESPEC{:},filled{:},varargin{:},'Color',CMAP(1,:),'Visible','off');
 end
else
 % do nothing
 if do2d
  quiver(AX,NaN,NaN,1,1,S,...
   LINESPEC{:},filled{:},varargin{:},'Visible','off'); 
 else
  quiver3(AX,NaN,NaN,NAN,1,1,1,S,...
   LINESPEC{:},filled{:},varargin{:},'Visible','off');
 end
end
% Bands between.
hold(AX,'on')
for k = 2:Nc-1
 ind = (L>=(Lmin+dL*(k-1))) & (L<(Lmin+dL*k));
 if ~isempty(ind)
  if do2d
   HQ(k) = quiver(AX,x(ind),y(ind),U(ind),V(ind),S,...
      LINESPEC{:},filled{:},varargin{:},'Color',CMAP(k,:),'Visible','off');
  else
   HQ(k) = quiver3(AX,x(ind),y(ind),z(ind),U(ind),V(ind),W(ind),S,...
      LINESPEC{:},filled{:},varargin{:},'Color',CMAP(k,:),'Visible','off');
  end
 end
end
% Last band
if Nc==1, ind = [];
else      ind = (L>=(Lmin+dL*(Nc-1))); end
if ~isempty(ind)
 if do2d
  HQ(Nc) = quiver(x(ind),y(ind),U(ind),V(ind),S,...
     LINESPEC{:},filled{:},varargin{:},'Color',CMAP(Nc,:),'Visible','off');
 else
  HQ(Nc) = quiver3(x(ind),y(ind),z(ind),U(ind),V(ind),W(ind),S,...
     LINESPEC{:},filled{:},varargin{:},'Color',CMAP(Nc,:),'Visible','off');
 end
end

% Finishes to draw arrows.  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Good handles.
igood = ishandle(HQ);
HQ    = HQ(igood);
if isempty(HQ)
 if nargout<1, clear(HQ)
 else HCB = []; HL  = []; end
 return
end
CMAP = CMAP(igood,:);
Nc   = size(CMAP,1);

% Draws lines for 'feathers'.
if ~strcmp('quiver',TYPE)
 % Looks for '-d'
 ind = strfind(lower(TYPE),'-d');
 if (length(TYPE)>1) && ~isempty(ind)
  dodate = true; TYPE = TYPE(1:ind-1);
 else
  dodate = false;
 end
 if strncmpi('feather',TYPE,length(TYPE)), FLINESPEC = {};
 else                                      FLINESPEC = {TYPE}; end
 try
  if do2d
    HL = plot([NaN(su(1),1), x, NaN(su(1),1)].',...
              [NaN(su(1),1), y, NaN(su(1),1)].',...
              FLINESPEC{:},'Parent',AX);
             % If an error occurs maybe is because 
             % TYPE was not a valid LINESPEC.
  else
   HL = plot3([NaN(su(1),1), x, NaN(su(1),1)].',...
              [NaN(su(1),1), y, NaN(su(1),1)].',...
              [NaN(su(1),1), z, NaN(su(1),1)].',...
              FLINESPEC{:},'Parent',AX);
            % If an error occurs maybe is because 
            % TYPE was not a valid LINESPEC.
  end
 catch
  error('CVARGAS:quiver2:unrecognizedLinespecFeather',...
   'Unrecognized ''t='',LINESPEC for the feather horizontal line.')
 end
else  HL = [];
end

if ~ihold, hold(AX,'off'), end

% Widens line.
if Nc==1, set(HQ(1),'LineWidth',WLIM(1))
else
 for k = 1:Nc, set(HQ(k),'LineWidth',interp1([1 Nc],WLIM,k)), end
end

% Changes colormap and limits.
if Nc>1
 if autoclim, caxis(AX,[Lmin Lmax]), end
 colormap(AX,CMAP) % forces COLORBAR generation
 HCB = colorbar('Peer',AX);
else HCB = [];
end

% Changes SLANT.
if strcmp('screen',SLANT), axis(AX,'equal'), end

% Changes Limits.
axis(AX,'tight'), xyz = axis(AX);
if ~autoxlim, set(AX,'XLim' ,limx ) 
else
 dx = xyz(2)-xyz(1); set(AX,'XLim',xyz(1:2)+K*dx*[-1 1])
 if autoxtic && strcmp('quiver',TYPE)
  set(AX,'XTickMode','auto','XTickLabelMode','auto')
 end
end
if ~autoylim, set(AX,'YLim' ,limy )
else
 dy = xyz(4)-xyz(3); set(AX,'YLim',xyz(3:4)+K*dy*[-1 1])
 if autoytic, set(AX,'YTickMode','auto','YTickLabelMode','auto'), end
end
if ~autozlim, set(AX,'ZLim' ,limz )
elseif ~do2d && (length(xyz)==6)
 dz = xyz(6)-xyz(5); set(AX,'XLim',xyz(5:6)+K*dz*[-1 1])
 if autoztic, set(AX,'ZTickMode','auto','ZTickLabelMode','auto'), end
end
if ~autoxtic, set(AX,'XTick',xtick), end
%if ~autoytic, set(AX,'YTick',ytick), end
%if ~autoztic, set(AX,'ZTick',ztick), end

% Sets feather x-ticks.
if ~strcmp(TYPE,'quiver') && ~all(x(:)==(1:numel(x))')
 if autoxlim
  limx = get(AX,'XLim');
  if all(x(~inan(:))==1), limx2  = (limx-1)+xlim(1);
  else limx2 = interp1([1 su(2)],xlim,limx,'linear','extrap'); end
 end
 if autoxtic
  xtick = get(AX,'XTick');
  if all(x(~inan(:))==1), xtick2 = (xtick-1)+xlim(1);
  else xtick2 = interp1([1 su(2)],xlim,xtick,'linear','extrap'); end
 end
 tempa = gca;
 ha2 = axes(...
  'HandleVisibility','off',...
  'Visible'         ,'off',...
  'Parent'          ,ancestor(AX,{'figure','uipanel'}),...
  'Units'           ,get(AX,'Units'),...
  'Position'        ,get(AX,'Position'),...
  'XLim'            ,limx2);
 if ~autoxtic, set(ha2,'XTick',xtick2), end
 if dodate
  if autoxtic
   if exist('tlabel.m','file')==2, tlabel(ha2,'x','keeplimits')
   else datetick(ha2,'x','keeplimits'), end
  else
   if exist('tlabel.m','file')==2, tlabel(ha2,'x','keeplimits','keepticks')
   else datetick(ha2,'x','keeplimits','keepticks'), end
  end
 end
 if autoxtic
  if all(x(~inan(:))==1), xtick = get(ha2,'XTick')-xlim(1)+1;
  else
   xtick = interp1(xlim,[1 su(2)],get(ha2,'XTick'),'linear','extrap');
  end
 end
 set(AX,'XTick',xtick,'XTickLabel',get(ha2,'XTickLabel'))
 % In case TLABEL was used.
 temp = get(get(ha2,'XLabel'),'String');
 if ~isempty(temp) && isempty(get(get(AX,'XLabel'),'String'))
  set(get(AX,'XLabel'),'String',temp);
 end
 delete(ha2)
 axes(tempa)
end

% Do ARROWs.
if doarrows
 HQ = mat2cell(HQ,ones(Nc,1),1);
 % Gets LINESPEC properties.
 if ~isempty(LINESPEC)
  if ~isempty(LineStyle), LineStyle = {'LineStyle',LineStyle};
  else                    LineStyle = {}; end
  if ~isempty(Marker),    Marker = {'Marker',Marker};
  else                    Marker = {};    end
 else
  LineStyle = {}; Marker = {};
 end
 % It need to create an invisible axes because ARROW functions does not
 % accept 'Parent' property.
 tempf = get(0    ,'CurrentFigure'); % Gets current figure handle
 tempa = get(tempf,'CurrentAxes');   % Gets current axes handle
 AF    = ancestor(AX,{'figure','uipanel'});
 hf2 = figure(...                    % Generates invisible figure
  'Visible' ,'on',...
  'Units'   ,get(AF,'Units'),...
  'Position',get(AF,'Position'));
 ha2 = axes(...                      % Generates invisible axes
  'Parent'  ,hf2,...
  'Units'   ,get(AX,'Units'),...
  'Position',get(AX,'Position'),...
  'XLim'    ,get(AX,'XLim'),...
  'YLim'    ,get(AX,'YLim'),...
  'ZLim'    ,get(AX,'ZLim'));
 set(ha2,'Units','Pixels')
 l2p = get(ha2,'Position');
 l2p = min(l2p(3:4)./[diff(get(ha2,'XLim')) diff(get(ha2,'YLim'))]);
 set(ha2,'Units',get(AX,'Units'))
 hold(ha2,'on')
 kbad = false(Nc,1);
 lw   = zeros(Nc,1);
 if strcmpi('Edge',EorF), ForE = 'Face'; else ForE = 'Edge'; end
 for k = 1:Nc
  % Start to draw invisible ARROWs
  HL = get(HQ{k},'Children');
  xa = get(HL(1),'XData'); xa = xa(:);
  ya = get(HL(1),'YData'); ya = ya(:);
  za = get(HL(1),'ZData'); za = za(:);
  lw(k) = get(HL(1),'LineWidth');
  delete(HQ{k})
  if isempty(xa) || isempty(ya)
   kbad(k) = true;
   continue
  end
  if isempty(za)
   Start = [xa(1:3:end) ya(1:3:end)];
   Stop  = [xa(2:3:end) ya(2:3:end)];
   u     = hypot(Stop(:,1)-Start(:,1),Stop(:,2)-Start(:,2));
   u     = max(u);
  else
   Start = [xa(1:3:end) ya(1:3:end) za(1:3:end)];
   Stop  = [xa(2:3:end) ya(2:3:end) za(2:3:end)];
   u     = hypot(Stop(:,1)-Start(:,1),Stop(:,2)-Start(:,2));
   u     = hypot(u,Stop(:,3)-Start(:,3));
   u     = max(u);
  end
  if u<eps*10^6
   kbad(k) = true;
   continue
  end
  if isempty(filled)
   HQ{k} = arrow(Start,Stop,'Length',alpha*l2p*u,Aopt{:},...
    [EorF 'Color'],CMAP(k,:),[ForE 'Color'],'none',...
    LineStyle{:}, Marker{:},varargin{:});
  else
   HQ{k} = arrow(Start,Stop,'Length',alpha*l2p*u,Aopt{:},...
    'FaceColor',CMAP(k,:),'EdgeColor',CMAP(k,:),...
    LineStyle{:}, Marker{:},varargin{:});
  end
 end
 clear xa ya za Start Stop
 lw(kbad) = [];
 HQ(kbad) = [];
 Nc       = length(HQ);
 % Draws the ARROWS on AX:
 ihold = ishold(AX);
 hold(AX,'on')
 for k = 1:Nc
  set(HQ{k},...
   {'Parent'}   ,{AX},...
   {'LineWidth'},{lw(k)},...
   {'Clipping'},{'on'})
 end
 if ~ihold, hold(AX,'off'), end
 delete(hf2)
 figure(tempf)
 axes(tempa)
else
 % USE QUIVER instead of ARROW.
 for k = 1:Nc, set(HQ(k),'Visible',visible), end
end


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Check nargout.
if nargout==0, clear HQ, end


% [EOF]   quiver2.m