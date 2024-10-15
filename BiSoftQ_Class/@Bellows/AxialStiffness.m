function k_bell = AxialStiffness(obj,stroke)

h = (obj.L+stroke)/(2*obj.ns);
alpha_s = asin(h/obj.AB);
lhs = h/tan(alpha_s);
D = obj.d+2*lhs;

k_bell = obj.E/(2*(1-obj.nu^2))*pi*(D+obj.d)/2*obj.t^3/lhs^3*1/obj.ns;

end