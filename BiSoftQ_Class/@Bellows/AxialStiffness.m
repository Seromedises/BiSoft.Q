function k_bell = AxialStiffness(obj,stroke)

h = (obj.L+stroke)/(2*obj.nb);
alpha_b = asin(h/obj.AB);
lhb = h/tan(alpha_b);
D = obj.d+2*lhb;

k_bell = obj.E/(2*(1-obj.nu^2))*pi*(D+obj.d)/2*obj.t^3/lhb^3*1/obj.nb;

end