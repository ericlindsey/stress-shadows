

N=128;
dx=5;
x=(-N/2:N/2-1)*dx;
y=(-N/2:N/2-1)*dx;

[xx,yy]=meshgrid(x,y);


depth=80;
L=200;
W=100;
[ux,uy,uz]=unicycle.greens.computeOkada85(1,xx(:),yy(:),0.25,pi/2,depth,L,W,1,0);

down=4;
index=reshape((1:numel(xx)),size(xx));
index=index(1:down:end,1:down:end);
index=index(:);


figure(110);clf;
subplot(2,2,1);cla;
scatter(xx(:),yy(:),[],ux(:),'filled');
colorbar
axis equal tight, box on
subplot(2,2,2);cla;
scatter(xx(:),yy(:),[],uy(:),'filled');
colorbar
axis equal tight, box on
subplot(2,2,3);cla;
scatter(xx(:),yy(:),[],uz(:),'filled');
colorbar
axis equal tight, box on
subplot(2,2,4);cla;hold on
scatter(xx(:),yy(:),[],uz(:),'filled');
quiver(xx(index),yy(index),ux(index)',uy(index)')
colorbar
axis equal tight, box on