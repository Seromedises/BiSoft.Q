function [X,Y,Z,Nx,Ny,Nz] = BellowsSurface(obj,alpha_b,h,rcvb,Rvb,Rt1,h1,rccb,Rcb,Rt2,stroke)

L = obj.L;
nb = obj.nb;
section_view = obj.section_view;
angle1_section_view = obj.angle1_section_view;
angle2_section_view = obj.angle2_section_view;

A1A2 = (h-Rcb*cos(alpha_b)-h1)/sin(alpha_b);
rho_v_I = sqrt((rcvb^2+Rvb^2-2*rcvb*Rvb*sin(alpha_b)));

v_I = asin((Rvb*cos(alpha_b))/rho_v_I);
v_II = asin((h-Rcb*cos(alpha_b))/sqrt(rho_v_I^2+A1A2^2-2*rho_v_I*A1A2*cos(pi-(alpha_b-v_I))));
v_III = atan(h/(rccb+Rcb));

rho_func1 = @(v) rcvb*cos(v)-sqrt(Rvb^2-(rcvb*sin(v)).^2);
rho_func2 = @(v) (Rvb*cos(alpha_b)-Rt1*tan(alpha_b))./(sin(v)-tan(alpha_b)*cos(v));
rho_func3 = @(v) rccb/(sin(atan(rccb/h)))*sin(v+atan(rccb/h))+sqrt(Rcb^2-(rccb/(sin(atan(rccb/h)))*cos(v+atan(rccb/h))).^2);

rho_func = @(v) rho_func1(v).*(0<=v & v<=v_I)+...
                rho_func2(v).*(v_I<v & v<=v_II)+...
                rho_func3(v).*(v_II<v & v<=v_III);

r_func = @(u,v) rho_func(v).*cos(v);
z_func = @(u,v) rho_func(v).*sin(v)-(L/2+stroke)*ones(size(v));

discr_u = 40; 
discr_v1 = 3;
discr_v2 = 2;
discr_v3 = 3;

u = linspace(0,2*pi,discr_u)';
v1 = linspace(0,v_I,discr_v1)';
v2 = linspace(v_I,v_II,discr_v2)';
v3 = linspace(v_II,v_III,discr_v3)';
v = [v1; v2; v3];

[U,V] = meshgrid(u,v);

R = r_func(U,V);
Z = z_func(U,V);

[R_mirror,Z_mirror] = MirrorPoints(obj,0,R,Z-h+(L/2+stroke));

R_tmp = [R; R_mirror];
Z_tmp = [Z; Z_mirror+h-(L/2+stroke)];

R = NaN(size(R_tmp,1)*nb,size(R_tmp,2));
Z = NaN(size(Z_tmp,1)*nb,size(Z_tmp,2));

R(1:size(R_tmp,1),:) = R_tmp;
Z(1:size(Z_tmp,1),:) = Z_tmp;

for ii=2:nb
    R(((ii-1)*size(R_tmp,1)+1):size(R_tmp,1)*ii,:) = R_tmp;
    Z(((ii-1)*size(Z_tmp,1)+1):size(Z_tmp,1)*ii,:) = Z_tmp+(ii-1)*2*h;
end

theta = repmat(U,2*nb,1);

if section_view==1
    [~,col_indx_angle1] = find(theta(1,:)<=angle1_section_view,1,'last');
    [~,col_indx_angle2] = find(theta(1,:)<=angle2_section_view,1,'last');

    if isempty(col_indx_angle1)
        col_indx_angle1 = 1;
    end

    if isempty(col_indx_angle2)
        col_indx_angle2 = numel(u);
    end

    theta = theta(:,col_indx_angle1:col_indx_angle2);
    R = R(:,col_indx_angle1:col_indx_angle2);
    Z = Z(:,col_indx_angle1:col_indx_angle2);

end

X = R.*cos(theta);
Y = R.*sin(theta);


Nx1 = @(u,v) -cos(u).*cos(v).*(sin(v).*(rcvb.*sin(v)-rcvb.^2.*cos(v).*sin(v).*1.0./sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))-cos(v).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2));
Nx2 = @(u,v) (cos(u).*cos(v).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)).*((cos(v).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)))./(sin(v)-tan(alpha_b).*cos(v))-sin(v).*(cos(v)+tan(alpha_b).*sin(v)).*1.0./(sin(v)-tan(alpha_b).*cos(v)).^2.*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b))))./(sin(v)-tan(alpha_b).*cos(v));
Nx3 = @(u,v) cos(u).*cos(v).*(sin(v).*(h.*cos(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)+h.^2.*cos(v+atan(rccb./h)).*sin(v+atan(rccb./h)).*1.0./sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0)).*(1.0./h.^2.*rccb.^2+1.0))+cos(v).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0))).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0));

Ny1 = @(u,v) -cos(v).*sin(u).*(sin(v).*(rcvb.*sin(v)-rcvb.^2.*cos(v).*sin(v).*1.0./sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))-cos(v).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2));
Ny2 = @(u,v) (cos(v).*sin(u).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)).*((cos(v).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)))./(sin(v)-tan(alpha_b).*cos(v))-sin(v).*(cos(v)+tan(alpha_b).*sin(v)).*1.0./(sin(v)-tan(alpha_b).*cos(v)).^2.*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b))))./(sin(v)-tan(alpha_b).*cos(v));
Ny3 = @(u,v) cos(v).*sin(u).*(sin(v).*(h.*cos(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)+h.^2.*cos(v+atan(rccb./h)).*sin(v+atan(rccb./h)).*1.0./sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0)).*(1.0./h.^2.*rccb.^2+1.0))+cos(v).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0))).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0));

Nz1 = @(u,v) cos(u).*cos(v).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2)).*(cos(u).*cos(v).*(rcvb.*sin(v)-rcvb.^2.*cos(v).*sin(v).*1.0./sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))+cos(u).*sin(v).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2)))+cos(v).*sin(u).*(cos(v).*sin(u).*(rcvb.*sin(v)-rcvb.^2.*cos(v).*sin(v).*1.0./sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))+sin(u).*sin(v).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2))).*(rcvb.*cos(v)-sqrt(-rcvb.^2.*sin(v).^2+Rvb.^2));
Nz2 = @(u,v) (cos(u).*cos(v).*((cos(u).*sin(v).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)))./(sin(v)-tan(alpha_b).*cos(v))+cos(u).*cos(v).*(cos(v)+tan(alpha_b).*sin(v)).*1.0./(sin(v)-tan(alpha_b).*cos(v)).^2.*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b))).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)))./(sin(v)-tan(alpha_b).*cos(v))+(cos(v).*sin(u).*((sin(u).*sin(v).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)))./(sin(v)-tan(alpha_b).*cos(v))+cos(v).*sin(u).*(cos(v)+tan(alpha_b).*sin(v)).*1.0./(sin(v)-tan(alpha_b).*cos(v)).^2.*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b))).*(Rvb.*cos(alpha_b)-Rt1.*tan(alpha_b)))./(sin(v)-tan(alpha_b).*cos(v));
Nz3 = @(u,v) -cos(u).*cos(v).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)).*(cos(u).*cos(v).*(h.*cos(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)+h.^2.*cos(v+atan(rccb./h)).*sin(v+atan(rccb./h)).*1.0./sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0)).*(1.0./h.^2.*rccb.^2+1.0))-cos(u).*sin(v).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)))-cos(v).*sin(u).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)).*(cos(v).*sin(u).*(h.*cos(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)+h.^2.*cos(v+atan(rccb./h)).*sin(v+atan(rccb./h)).*1.0./sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0)).*(1.0./h.^2.*rccb.^2+1.0))-sin(u).*sin(v).*(sqrt(Rcb.^2-h.^2.*cos(v+atan(rccb./h)).^2.*(1.0./h.^2.*rccb.^2+1.0))+h.*sin(v+atan(rccb./h)).*sqrt(1.0./h.^2.*rccb.^2+1.0)));

Nx_func = @(u,v) Nx1(u,v).*(0<=v & v<=v_I)+...
                 Nx2(u,v).*(v_I<v & v<=v_II)+...
                 Nx3(u,v).*(v_II<v & v<=v_III);

Ny_func = @(u,v) Ny1(u,v).*(0<=v & v<=v_I)+...
                 Ny2(u,v).*(v_I<v & v<=v_II)+...
                 Ny3(u,v).*(v_II<v & v<=v_III);

Nz_func = @(u,v) Nz1(u,v).*(0<=v & v<=v_I)+...
                 Nz2(u,v).*(v_I<v & v<=v_II)+...
                 Nz3(u,v).*(v_II<v & v<=v_III);

Nx = Nx_func(U,V);
Ny = Ny_func(U,V);
Nz = Nz_func(U,V);

Nx = [Nx; Nx];
Ny = [Ny; Ny];
Nz = [Nz; -Nz];

Nx = repmat(Nx,nb,1);
Ny = repmat(Ny,nb,1);
Nz = repmat(Nz,nb,1);















