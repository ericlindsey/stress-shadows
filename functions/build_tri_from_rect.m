function flt = build_tri_from_rect(rcv,params)
% function to convert rectangular fault mesh to triangles for cartesian coordinates
% INPUT
% src - unicycle rectangle source object
% params - structure containing the following:
% n1 - number of points to grid top and bottom of fault (generally set ~100)
% n2 - number of points to grid edge of fault (generally set ~10)
% delshift - distance to shift the top/bottom edge of the fault to allow for clean meshing (this is annoyingly manual, generally set as 2e3 (m))
% name - provide name of tri file output with no extension (palu_tri_xyz)
% trilen - target length for each triangular patch (m)
% OUTPUT
% flt - unicycle triangle source object
% Rishav Mallick, 2018, Earth Observatory of Singapore

import unicycle.*
% modify this so that you have mesh2d
addpath(genpath('/Users/elindsey/Dropbox/projects/stress_constraint_backslip/interseismic_stressinv/jointinv2/functions/mesh2d-master'))
%% set defaults
if exist('params')==0
    params.n1 = 100;
    params.n2 = 10;
    params.delshift = 0;
    params.name = 'dummy_tri_xyz';
    params.trilen = 2e3;
else
    if isfield(params,'n1')==0
        params.n1 = 100;
    end
    if isfield(params,'n2')==0
        params.n2 = 10;
    end 
    if isfield(params,'delshift')==0
        params.delshift = 0;
    end
    if isfield(params,'name')==0
        params.name = 'dummy_tri_xyz';
    end
    if isfield(params,'trilen')==0
        params.trilen = 2e3;
    end
end
        
%% locate top and bottom of rectangular fault mesh

Itop = rcv.x(:,3)==0;
Ibot = find_bottom(rcv);

bb_top_NZ = rcv.x(Itop,[2 3]);
bb_top_E = rcv.x(Itop,1); % for meshing we will use only 2 and 3, like in a projection -> view(90,0)
addW = [rcv.W(Ibot) rcv.W(Ibot) rcv.W(Ibot)].*rcv.dv(Ibot,:);
bb_bot_NZ = rcv.x(Ibot,[2 3]) - addW(:,[2 3]);
bb_bot_E = rcv.x(Ibot,1) - addW(:,1);

bb_bot_NZ=flipud(bb_bot_NZ);
bb_bot_E=flipud(bb_bot_E);

n1=params.n1;
n2=params.n2;

% add delshift to help spline use points that are at the same x or y coordinate
if bb_top_NZ(end,1)==bb_bot_NZ(1,1) || bb_bot_NZ(end,1)==bb_top_NZ(1,1)
    delshift = params.delshift;
else
    delshift = 0;
end
    
Y = [linspace(bb_top_NZ(1,1),bb_top_NZ(end,1),n1) ...
    linspace(bb_top_NZ(end,1),bb_bot_NZ(1,1)-delshift,n2) ...
    linspace(bb_bot_NZ(1,1),bb_bot_NZ(end,1),n1) ...
    linspace(bb_bot_NZ(end,1)+delshift,bb_top_NZ(1,1),n2)]';

Z = [spline(bb_top_NZ(1:end,1),bb_top_NZ(1:end,2),Y(1:n1));... 
    spline([bb_top_NZ(end,1) bb_bot_NZ(1,1)-delshift],[bb_top_NZ(end,2) bb_bot_NZ(1,2)],Y(n1+1:n1+n2));...
    spline(bb_bot_NZ(1:end,1),bb_bot_NZ(1:end,2),Y(n1+n2+1:2*n1+n2));...
    spline([bb_bot_NZ(end,1)+delshift bb_top_NZ(1,1)],[bb_bot_NZ(end,2) bb_top_NZ(1,2)],Y(2*n1+n2+1:end))];

%% do actual meshing
X = [interp1(bb_top_NZ(:,1),bb_top_E,Y(1:n1),'cubic');...
    linspace(bb_top_E(end),bb_bot_E(1),n2)';...
    interp1(bb_bot_NZ(:,1),bb_bot_E,Y(n1+n2+1:2*n1+n2),'cubic');...
    linspace(bb_bot_E(end),bb_top_E(1),n2)'];

p = [Z Y];
edge = [];
for i = 1:length(Z)
    if i ~= length(Z)
        edge(i,1:2) = [i i+1];
    else
        edge(i,1:2) = [i 1];
    end
end
% build lowest resolution possible
% [vert,etri,tria,tnum] = refine2(p,edge) ;

hfun = params.trilen ;            % uniform "target" edge-lengths
[vert,etri,tria,tnum] = refine2(p,edge,[],[],hfun) ;

% here the vertices are ordered in terms of N(y) and Z(z) locations. Need
% to switch this around for the actual fault mesh

% provide the values of known (x,y,z) points ON the fault. The most
% convenient way would be to simply use the provided slab geometry from
% slab1.0 but it needs to shifted to the free surface. Instead, what I do
% is pin the trace to 0, and then use interpolation to the down-dip extent
% at 50-60 km depth. For greater accuracy, I can add in the 20km and 40km
% depth contours as well to (xf,yf,zf)
zf = [Z;rcv.xc(:,3)];
yf = [Y;rcv.xc(:,2)];
xf = [X;rcv.xc(:,1)];

xint = griddata(yf,zf,xf,vert(:,2),vert(:,1),'cubic');
TF = isnan(xint);% | xint < 0;
xint(TF) = 0;
% vert(:,1) - depth
% vert(:,2) - Y(north)
% xint - X(east)
%% write flt file

filenameout = params.name;
write_file = 1;
if write_file ==1
    nedfname=[filenameout '.ned'];
    fileID = fopen(nedfname,'w');
    fprintf(fileID,'%s\n',        '#n  X1   X2    X3   \n');
    for i = 1:length(vert(:,1))
        fprintf(fileID,'%.0f %.4f %.4f %.4f \n',...
            i, vert(i,2), xint(i), -vert(i,1));
    end
    fclose(fileID);
    
    trifname=[filenameout '.tri'];
    fileID = fopen(trifname,'w');
    fprintf(fileID,'%s\n',        '#n  Id 1   Id 2    Id 3     rake(deg)');
    for i = 1:length(tria(:,1))
        fprintf(fileID,'%4.0f %4.0f %4.0f %4.0f %4.2f \n',...
            i, tria(i,2), tria(i,1), tria(i,3), 90);
    end
    fclose(fileID);
end

G=30e9;
nu=1/4;

flt=geometry.triangleReceiver(filenameout,greens.nikkhoo15(G,nu));

end