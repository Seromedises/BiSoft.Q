function F_pull_approx = PullForceApprox(obj,GeometricParams,p,stroke,method,approx_method)

L = obj.L;
N = obj.N;
alpha = obj.alpha;
Rie = obj.Rie;
Roe = obj.Roe;
[~,rc,phi_c,phi_v,~,~,~,~,~] = struct('x', num2cell(GeometricParams)).x;

rv = obj.lv/(2*phi_v);
Rom = obj.Roe+rc*(1-cos(phi_c));
Rim = obj.Rie+rv*(1-cos(phi_v));

if strcmp(method,'method_1')
    phi_angle = phi_v;
else
    phi_angle = phi_c;
end

if strcmp(approx_method,'approx_1')

    F_pull_approx = p*N*(sin(alpha/2)*(L+stroke)*Rie)/tan(phi_angle);

elseif strcmp(approx_method,'approx_2')

    b = sqrt((Rie*sin(alpha/2))^2+(Roe-Rie*cos(alpha/2))^2);
    B = sqrt((Rim*sin(alpha/2))^2+(Rom-Rim*cos(alpha/2))^2);
    epsilon = atan((Roe-Rie*cos(alpha/2))/(Rie*sin(alpha/2)));

    if strcmp(method,'method_1')
        angle = epsilon-alpha/2;
    else
        angle = epsilon;
    end

    Ai = 1/4*(b+B)*(L+stroke); % area of the projected trapezium

    Fp_radial = 2*(2*p*Ai*cos(angle));
    F_pull_approx = N*Fp_radial/(2*tan(phi_angle));

end


end
