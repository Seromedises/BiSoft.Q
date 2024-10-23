function volume = Volume(obj,GeometricParams)

nb = obj.nb;
[alpha_b,h,rcvb,Rvb,Rt1,h1,rccb,Rcb,Rt2,h2] = struct('x', num2cell(GeometricParams)).x;

Ab1 = 0.5*(Rvb^2*(pi/2-alpha_b)-h1*Rvb*cos(pi/2-alpha_b));
Ab2 = 0.5*(Rcb^2*(pi/2-alpha_b)-(h-h2)*Rcb*cos(pi/2-alpha_b));

V1 = pi*Rt1^2*h1-2*pi*rcvb*Ab1;
V2 = 1/3*pi*(h2-h1)*(Rt1^2+Rt1*Rt2+Rt2^2);
V3 = pi*Rt2^2*(h-h2)+2*pi*rccb*Ab2;

volume = (V1+V2+V3)*2*nb;

end