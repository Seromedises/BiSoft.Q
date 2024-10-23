function PlotNominalGeometry(obj)

colors = [0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980];

% initial configuration
[~,~,~,~,~,X_profile,Y_profile,Z_profile,~,~,~] = PleatedMembraneSurface(obj,obj.rc0,obj.rc,obj.rs1,obj.x0s1,obj.rs2,obj.x0s2,obj.vmax);

% 3d plot of the parametric surface
fig(1) = figure('Name','Surface');
surf(X_profile,Y_profile,Z_profile)
shading interp
hold on
theta = linspace(0,2*pi,100);
plot3(obj.Rib*cos(theta),obj.Rib*sin(theta),Z_profile(1,1)*ones(size(theta)),'k-.')
plot3(obj.Rib*cos(theta),obj.Rib*sin(theta),Z_profile(1,25)*ones(size(theta)),'k-.')
plot3(obj.Rib*cos(theta),obj.Rib*sin(theta),Z_profile(1,end)*ones(size(theta)),'k-.')
hold off
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
view(35,30)
title("Pleated membrane")

% plot the complete profile on middle cross section
middle_idx = size(X_profile,2)/2;
fig(2) = figure('Name','Middle and end plate cross sections');
hold on
plot(X_profile(:,middle_idx),Y_profile(:,middle_idx),'b')
plot(X_profile(:,end),Y_profile(:,end),'r')
theta = linspace(0,2*pi,100);
plot(obj.Rib*cos(theta),obj.Rib*sin(theta),'k-.')
legend('membr. (middle)','membr. (end plate)','bellows','location','southwest')
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
box on
hold off
title("Cross sections")

% plot of the longitudinal cross sections
fig(3) = figure('Name','Longitudinal cross sections');
theta_v = obj.phi_v*linspace(-1,1,100);
xv = obj.rv*cos(theta_v)+obj.Rim-obj.rv;
yv = obj.rv*sin(theta_v);
theta_c = obj.phi_c*linspace(-1,1,100);
xc = obj.rc*cos(theta_c)+obj.Rom-obj.rc;
yc = obj.rc*sin(theta_c);
hold on
P3 = plot(obj.Rib*[1 1 -1 -1 1],obj.L/2*[-1 1 1 -1 -1],'k-.');
P1 = plot([-xc(1) xc -xc(end) -flip(xc)],[yc(1) yc yc(end) flip(yc)],'Color',colors(1,:));
P2 = plot([-xv(1) xv -xv(end) -flip(xv)],[yv(1) yv yv(end) flip(yv)],'Color',colors(2,:));
axis equal
hold off
box on
xlim([-30 30])
ylim([-25 25])
xlabel('(mm)')
ylabel('Z (mm)')
legend([P1 P2],'crest','valley')
title('Longitudinal cross sections')

end