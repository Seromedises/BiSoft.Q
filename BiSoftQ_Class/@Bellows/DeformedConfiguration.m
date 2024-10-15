function DeformedConfigParameters = DeformedConfiguration(obj,stroke)

h = (obj.L+stroke)/(2*obj.ns);
alpha_s = asin(h/obj.AB);
lhs = h/tan(alpha_s);
Rs1 = (obj.Rib-obj.d/2)*sin(alpha_s)/(1-sin(alpha_s));
h1 = Rs1*cos(alpha_s);
Rs2 = (obj.ltot-Rs1*(pi/2-alpha_s)-(h-h1)/sin(alpha_s))/(pi/2-alpha_s-1/tan(alpha_s));
h2 = h-Rs2*cos(alpha_s);
AA1 = Rs1/tan(alpha_s);
A2B = Rs2/tan(alpha_s);
Rt1 = obj.d/2+AA1*cos(alpha_s);
Rt2 = obj.d/2+(obj.AB-A2B)*cos(alpha_s);
rcs1 = obj.Rib+Rs1;
rcs2 = obj.d/2+lhs-Rs2/sin(alpha_s);

DeformedConfigParameters = [alpha_s,h,rcs1,Rs1,Rt1,h1,rcs2,Rs2,Rt2,h2];

end
