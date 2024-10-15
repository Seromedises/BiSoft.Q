function volume = Volume(obj,GeometricParams)

ns = obj.ns;
[alpha_s,h,rcs1,Rs1,Rt1,h1,rcs2,Rs2,Rt2,h2] = struct('x', num2cell(GeometricParams)).x;

Ab1 = 0.5*(Rs1^2*(pi/2-alpha_s)-h1*Rs1*cos(pi/2-alpha_s));
Ab2 = 0.5*(Rs2^2*(pi/2-alpha_s)-(h-h2)*Rs2*cos(pi/2-alpha_s));

V1 = pi*Rt1^2*h1-2*pi*rcs1*Ab1;
V2 = 1/3*pi*(h2-h1)*(Rt1^2+Rt1*Rt2+Rt2^2);
V3 = pi*Rt2^2*(h-h2)+2*pi*rcs2*Ab2;

volume = (V1+V2+V3)*2*ns;

end