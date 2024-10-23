function PlotNominalGeometry(obj)

[X,Y,Z] = BellowsSurface(obj,obj.alpha_b,obj.h,obj.rcvb,obj.Rvb,obj.Rt1,obj.h1, ...
    obj.rccb,obj.Rcb,obj.Rt2,0);

fig(1) = figure('Name','Bellows surface - initial configuration');
hold on
surf(X,Y,Z,'EdgeColor','none','FaceColor',[0.90 0.90 0.90],'FaceAlpha',0.9)
plot3(X(1,:),Y(1,:),Z(1,:),'k-','LineWidth',1.5)
for ii=2:2*obj.nb  
    plot3(X((ii-1)*size(Z,1)/(2*obj.nb),:),Y((ii-1)*size(Z,1)/(2*obj.nb),:),Z((ii-1)*size(Z,1)/(2*obj.nb),:),'k-','LineWidth',1.5)
end
plot3(X(end,:),Y(end,:),Z(end,:),'k-','LineWidth',1.5)
if obj.section_view==1
    plot3(X(:,1),Y(:,1),Z(:,1),'k-','LineWidth',1.5)
    plot3(X(:,end),Y(:,end),Z(:,end),'k-','LineWidth',1.5)
end
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
view(35,30)
grid off
hold off

fig(2) = figure('Name','Bellows longitudinal cross section');
R_sec = sqrt(X(:,1).^2+Y(:,1).^2);
[R_sec_mirror,Z_sec_mirror] = MirrorPoints(obj,pi/2,R_sec,Z(:,1));
R_sec = [R_sec;R_sec_mirror;R_sec(1)];
Z_sec = [Z(:,1);Z_sec_mirror;Z(1,1)];
plot(R_sec,Z_sec,'k-')
axis equal
max_R = max(R_sec,[],'all');
xlim(ceil(max_R/10)*10*[-1 1])
ylim([-obj.L/2-5 obj.L/2+5])
xlabel('X (mm)')
ylabel('Z (mm)')

end