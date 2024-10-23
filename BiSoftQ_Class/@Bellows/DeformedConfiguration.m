function DeformedConfigParameters = DeformedConfiguration(obj,stroke)

h = (obj.L+stroke)/(2*obj.nb);
alpha_b = asin(h/obj.AB);
lhs = h/tan(alpha_b);
Rvb = (obj.Rib-obj.d/2)*sin(alpha_b)/(1-sin(alpha_b));
h1 = Rvb*cos(alpha_b);
Rcb = (obj.ltot-Rvb*(pi/2-alpha_b)-(h-h1)/sin(alpha_b))/(pi/2-alpha_b-1/tan(alpha_b));
h2 = h-Rcb*cos(alpha_b);
AA1 = Rvb/tan(alpha_b);
A2B = Rcb/tan(alpha_b);
Rt1 = obj.d/2+AA1*cos(alpha_b);
Rt2 = obj.d/2+(obj.AB-A2B)*cos(alpha_b);
rcs1 = obj.Rib+Rvb;
rcs2 = obj.d/2+lhs-Rcb/sin(alpha_b);

DeformedConfigParameters = [alpha_b,h,rcs1,Rvb,Rt1,h1,rcs2,Rcb,Rt2,h2];

end
