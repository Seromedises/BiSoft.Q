function [X,Y,Z,X_curve,Y_curve,X_profile,Y_profile,Z_profile] = PleatedMembraneSurface(obj,rc0,rc,rs1,x0s1,rs2,x0s2,vmax)

N = obj.N;
alpha = obj.alpha;
delta = obj.delta;
gamma = obj.gamma;

% parametric surface
% define equations of the parametric surface
r = @(v) rc0*cos(v)+sqrt(rc^2-(rc0*sin(v)).^2);

os1 = @(v) sqrt(rs1^2-(r(v).*sin(v)).^2)+x0s1;
os2 = @(v) sqrt(rs2^2-(r(v).*sin(v)).^2)+x0s2;

beta_func = @(v) pi - atan2(os1(v)*cos(delta)-os2(v)*cos(alpha/2-gamma),os1(v)*sin(delta)-os2(v)*sin(alpha/2-gamma)); 

A_func = @(v) 1+sin(beta_func(v))/tan(delta)-cos(beta_func(v));
B_func = @(v) sin(beta_func(v)-alpha/2)/tan(gamma)+cos(beta_func(v)-alpha/2)-1;
C_func = @(v) (1-(B_func(v)+1).*cos(beta_func(v)-alpha/2))./(cos(beta_func(v)).*(1-A_func(v))-1); 

rho_func = @(u,v) (r(v).*cos(v).*((A_func(v)-1).*cos(u)+sqrt(1-((A_func(v)-1).*sin(u)).^2))./A_func(v)).*(0<=u & u<=delta) + ...
                  (r(v).*cos(v).*(1+(A_func(v)-1).*cos(beta_func(v)))./(A_func(v).*cos(u-beta_func(v)))).*(delta<u & u<=alpha/2-gamma) + ...
                  (r(v).*cos(v).*((B_func(v)+1).*cos(u-alpha/2)-sqrt(1-((B_func(v)+1).*sin(u-alpha/2)).^2))./(A_func(v).*C_func(v))).*(alpha/2-gamma<u & u<=alpha/2);

x_func = @(u,v) rho_func(u,v).*sin(u);
y_func = @(u,v) rho_func(u,v).*cos(u);
z_func = @(v) r(v).*sin(v);

% define intervals of the parameters u,v
discr_u = 25;
discr_v = 25;

v = linspace(0,vmax,discr_v)';
u = linspace(0,alpha/2,discr_u)';

[U,V] = meshgrid(u,v);

X = x_func(U,V)';
Y = y_func(U,V)';
Z = z_func(V)';

% mirror points
[X_mirror, Y_mirror] = MirrorPoints(obj,pi/2-alpha/2,X,Y);

X_curve = [X; X_mirror];
Y_curve = [Y; Y_mirror];

% repeat the curve to draw the complete profile
X_profile = zeros(size(X_curve,1)*(N+1),size(X_curve,2));
Y_profile = zeros(size(Y_curve,1)*(N+1),size(Y_curve,2));

for i=1:size(X_curve,2)
    [X_profile(:,i),Y_profile(:,i)] = RepeatProfile(obj,X_curve(:,i),Y_curve(:,i),N+1);
end

Z_profile = ones(size(X_profile,1),1)*Z(1,:);

X_profile = [flip(X_profile,2) X_profile];
Y_profile = [flip(Y_profile,2) Y_profile];
Z_profile = [-flip(Z_profile,2) Z_profile];

end
