function plot_gridded_slip_model(scenario)

N = scenario.sources{1}.geom.N;
slip=-scenario.modelVector(1:N);

% first, get the slip model values at the vertices, this will be
% non-uniqe
xout=zeros(N,1);
yout=zeros(N,1);
zout=zeros(N,1);
sout=zeros(N,1);

for itri = 1:N
    v = scenario.sources{1}.geom.vertices(itri);
    trislip = slip(itri);
    xtri = scenario.sources{1}.geom.x(scenario.sources{1}.geom.vertices(itri,:),1)/1e3;
    ytri = scenario.sources{1}.geom.x(scenario.sources{1}.geom.vertices(itri,:),2)/1e3;
    ztri = -scenario.sources{1}.geom.x(scenario.sources{1}.geom.vertices(itri,:),3)/1e3;
    xout=[xout; xtri; scenario.sources{1}.geom.xc(itri,1)/1e3];
    yout=[yout; ytri; scenario.sources{1}.geom.xc(itri,2)/1e3];
    zout=[zout; ztri; scenario.sources{1}.geom.xc(itri,3)/1e3];
    sout=[sout; repmat(trislip,4,1)];
end

% grid the x,y,z coordinates
yvec = -335:.5:358;
dvec=0:.1:20;


[yg,dg]=meshgrid(yvec,dvec);

% getting the x is harder
[ysort,Isort] = sort(scenario.sources{1}.geom.xc(:,2)/1e3);
[yuniq,Iuniq] = unique(ysort);
xsort=scenario.sources{1}.geom.xc(Isort,1)/1e3;
xuniq=xsort(Iuniq);
xvec = interp1(yuniq,xuniq, yvec,'nearest');
xg=repmat(xvec,size(yg,1),1);

[latvec,lonvec]=xy_to_latlon_polyconic(xg(:),yg(:),scenario.userParams.LOSlon0,scenario.userParams.LOSlat0);
lon_g=reshape(lonvec,size(xg));
lat_g=reshape(latvec,size(xg));

%

F = scatteredInterpolant(yout,zout,sout, 'linear');  % or 'linear', 'nearest'
% Now you can evaluate at any (xq, yq):
sg = F(yg, dg);

pcolor(lat_g,dg,sg),shading flat
colorbar
colormap(flipud(hot))
set(gca,'ydir','rev')

end