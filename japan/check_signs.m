% check all signs of the input files

%run('../functions/set_jointinv_path()')

scenario = Jointinv(0);
scenario.run_setup();

% flip sign because it was flipped earlier
%scenario.dataVector=-scenario.dataVector;

plateRake = scenario.sources{1}.rakeFile(:,1);
Vpl_ss = scenario.sources{1}.rakeFile(:,2);
Vpl_ds = scenario.sources{1}.rakeFile(:,3);


%%
figure(1),clf, hold on
scenario.sources{1}.geom.plotPatch(Vpl_ss)
scenario.sources{1}.geom.plotSlipVectors(Vpl_ss, 0*Vpl_ds,1000,'k')

pureRakeSlip=100*[...
    ones(scenario.sources{1}.geom.N, 1);...
    zeros(scenario.sources{1}.geom.N, 1)];

scenario.predVector = scenario.designMatrix * pureRakeSlip;

vecScale=1e3;
plot_coupling_vectors(scenario,vecScale)
daspect([ 1 1 1])



%% plot triangle cycle

figure(2),clf

% v = scenario.sources{1}.geom.vertices(1,:);
% vx=scenario.sources{1}.geom.x(v,:);
% % check the cross product
% v1=vx(2,:)-vx(1,:);
% v2=vx(3,:)-vx(2,:);
% cpcheck = cross(v1,v2);
% assert(cpcheck(3)<0,'Error: Vertices for all triangles in triangleReceiver object must be ordered in a clockwise sense.')

p1 = scenario.sources{1}.geom.x(scenario.sources{1}.geom.vertices(:,1),:);
p2 = scenario.sources{1}.geom.x(scenario.sources{1}.geom.vertices(:,2),:);
p3 = scenario.sources{1}.geom.x(scenario.sources{1}.geom.vertices(:,3),:);

v1 = p2 - p1;
v2 = p3 - p2;

cp=cross(v1,v2);

assert(max(cp(:,3))<=0,'Error: Vertices for all triangles in triangleReceiver object must be ordered in a clockwise sense.')


scatter(vx(:,1),vx(:,2),100*(1:3)',100*(1:3)','x')
caxis([0,400])
view(2)


