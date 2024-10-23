function [X,Y,Z] = BellowsSurface(obj,alpha_b,h,rcvb,Rvb,Rt1,h1,rccb,Rcb,Rt2,stroke)

L = obj.L;
nb = obj.nb;
section_view = obj.section_view;
angle1_section_view = obj.angle1_section_view;
angle2_section_view = obj.angle2_section_view;

theta = linspace(0,pi/2-alpha_b,24)';
r1 = rcvb-Rvb*cos(theta);
z1 = Rvb*sin(theta);
r2 = linspace(Rt1,Rt2,2)';
z2 = h1+(r2-Rt1)*tan(alpha_b);
r3 = rccb+Rcb*cos(flip(theta));
z3 = h-Rcb*sin(flip(theta));

r_vec = [r1;r2;r3];
z_vec = [z1;z2;z3];

[r_mirror,z_mirror] = MirrorPoints(obj,0,flip(r_vec),flip(z_vec)-h);

rs = [r_vec;flip(r_mirror)];
zs = [z_vec;flip(z_mirror)+h];

r = NaN(size(rs,1)*nb,1);
z = NaN(size(zs,1)*nb,1);

r(1:size(rs,1)) = rs;
z(1:size(zs,1)) = zs;

for ii=2:nb
    r(((ii-1)*size(rs,1)+1):size(rs,1)*ii) = rs;
    z(((ii-1)*size(rs,1)+1):size(rs,1)*ii) = zs+(ii-1)*2*h;
end

z = z-(L/2+stroke);

if section_view==0
    theta = linspace(0,2*pi,100);
else
    theta = linspace(angle1_section_view,angle2_section_view,100);
end
X = r.*cos(theta);
Y = r.*sin(theta);
Z = z.*ones(size(X));
