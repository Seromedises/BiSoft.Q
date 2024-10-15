function PlotNominalGeometry(obj)

L = obj.PleatedMembrane.L;
N = obj.PleatedMembrane.N;
Rib = obj.PleatedMembrane.Rib;
phi_v = obj.PleatedMembrane.phi_v;
rv = obj.PleatedMembrane.rv;
Rim = obj.PleatedMembrane.Rim;
phi_c = obj.PleatedMembrane.phi_c;
rc = obj.PleatedMembrane.rc;
Rom = obj.PleatedMembrane.Rom;
Roe = obj.PleatedMembrane.Roe;

section_view = obj.Bellows.section_view;

colors = [0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980];

% Pleated membrane
[~,~,~,~,~,X_membr,Y_membr,Z_membr] = obj.PleatedMembrane.PleatedMembraneSurface(...
    obj.PleatedMembrane.rc0,obj.PleatedMembrane.rc,obj.PleatedMembrane.rs1, ...
    obj.PleatedMembrane.x0s1,obj.PleatedMembrane.rs2,obj.PleatedMembrane.x0s2,obj.PleatedMembrane.vmax);

% Bellows
[X_bell,Y_bell,Z_bell] = obj.Bellows.BellowsSurface(obj.Bellows.alpha_s,obj.Bellows.h, ...
    obj.Bellows.rcs1,obj.Bellows.Rs1,obj.Bellows.Rt1,obj.Bellows.h1, ...
    obj.Bellows.rcs2,obj.Bellows.Rs2,obj.Bellows.Rt2,0);

if strcmp(obj.type,"A")
    FaceColorBell = [0.5 0.5 0.5];
    FaceAlphaBell = 1;
    FaceColorMembr = [0.90 0.90 0.90];
    FaceAlphaMembr = 0.9;
elseif strcmp(obj.type,"B")
    FaceColorBell = [0.90 0.90 0.90];
    FaceAlphaBell = 0.9;
    FaceColorMembr = [0.5 0.5 0.5];
    FaceAlphaMembr = 1;
end

% 3d plot of the parametric surface
fig(1) = figure('Name','Surface of the actuator (nominal configuration)');
% pleated membrane
surf(X_membr,Y_membr,Z_membr,'EdgeColor','none','FaceColor',FaceColorMembr,'FaceAlpha',FaceAlphaMembr)
hold on
for iii=1:2*N
    idx = (iii-1)*25+1;
    plot3(X_membr(idx,:),Y_membr(idx,:),Z_membr(idx,:),'k','LineWidth',1.5)
end
plot3(X_membr(:,1),Y_membr(:,1),Z_membr(:,1),'k','LineWidth',1.5)
% plot3(X_profile(:,25),Y_profile(:,25),Z_profile(:,25),'k','LineWidth',1.5)
plot3(X_membr(:,end),Y_membr(:,end),Z_membr(:,end),'k','LineWidth',1.5)
% bellows
surf(X_bell,Y_bell,Z_bell,'EdgeColor','none','FaceColor',FaceColorBell,'FaceAlpha',FaceAlphaBell)
plot3(X_bell(1,:),Y_bell(1,:),Z_bell(1,:),'k-','LineWidth',1.5)
for ii=2:size(Z_bell,1)/50   
    plot3(X_bell((ii-1)*50,:),Y_bell((ii-1)*50,:),Z_bell((ii-1)*50,:),'k-','LineWidth',1.5)
end
plot3(X_bell(end,:),Y_bell(end,:),Z_bell(end,:),'k-','LineWidth',1.5)
if section_view==1
    plot3(X_bell(:,1),Y_bell(:,1),Z_bell(:,1),'k-','LineWidth',1.5)
    plot3(X_bell(:,end),Y_bell(:,end),Z_bell(:,end),'k-','LineWidth',1.5)
end
hold off
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
view(35,30)
title(['BiSoft.Q type ',obj.type])

% plot the complete profile on middle cross section
middle_idx = size(X_membr,2)/2;
fig(2) = figure('Name','Middle and end plate cross sections');
hold on
plot(X_membr(:,middle_idx),Y_membr(:,middle_idx),'b')
plot(X_membr(:,end),Y_membr(:,end),'r')
theta = linspace(0,2*pi,100);
plot(Rib*cos(theta),Rib*sin(theta),'k-.')
legend('membr. (middle)','membr. (end plate)','bellows','interpreter','latex','location','southwest')
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
box on
hold off
title("Cross sections")

% plot of the longitudinal cross sections
fig(3) = figure('Name','Longitudinal cross sections');
% pleated membrane
theta_v = linspace(-phi_v,phi_v,100);
xv = rv*cos(theta_v)+Rim-rv;
yv = rv*sin(theta_v);
theta_c = linspace(-phi_c,phi_c,100);
xc = rc*cos(theta_c)+Rom-rc;
yc = rc*sin(theta_c);
hold on
P1 = plot([-xc(1) xc -xc(end) -flip(xc)],[yc(1) yc yc(end) flip(yc)],'Color',colors(1,:));
P2 = plot([-xv(1) xv -xv(end) -flip(xv)],[yv(1) yv yv(end) flip(yv)],'Color',colors(2,:));
% bellows
R_sec = sqrt(X_bell(:,1).^2+Y_bell(:,1).^2);
[R_sec_mirror,Z_sec_mirror] = obj.Bellows.MirrorPoints(pi/2,R_sec,Z_bell(:,1));
R_sec = [R_sec;R_sec_mirror;R_sec(1)];
Z_sec = [Z_bell(:,1);Z_sec_mirror;Z_bell(1,1)];
P3 = plot(R_sec,Z_sec,'k-');
if strcmp(obj.type,"A")
    n = 24;
    DD = Roe/n;
    plot([-Roe -Roe+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
    hold on
    for iii=2:n+1
        plot([-Roe+2*(iii-1)*DD -Roe+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
    end
elseif strcmp(obj.type,"B")
    n = 24;
    DD = Rib/n;
    plot([-Rib -Rib+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
    hold on
    for iii=2:n+1
        plot([-Rib+2*(iii-1)*DD -Rib+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
    end
end
axis equal
hold off
box on
max_R = max(R_sec,[],'all');
if strcmp(obj.type,"A")
    xlim(ceil(Rom/10)*10*[-1 1])
elseif strcmp(obj.type,"B")
    xlim(ceil(max_R/10)*10*[-1 1])
end
ylim([-L/2-5 L/2+5])
xlabel('(mm)')
ylabel('Z (mm)')
legend([P1 P2 P3],'crest','valley','bellows','interpreter','latex')
title('Longitudinal cross sections')

end

