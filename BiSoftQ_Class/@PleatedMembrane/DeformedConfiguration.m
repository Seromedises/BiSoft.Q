function [DeformedConfigParameters,sol,StopFlag] = DeformedConfiguration(obj,stroke,initial_guess)

StopFlag = false;

% find angle and radius of the crest arc in the deformed configuration 
phi_c = fzero(@(phi) sincc(phi/pi)-(obj.L+stroke)/obj.lc,[0 pi]);
rc = obj.lc/(2*phi_c);

% find angle and radius of the valley arc in the deformed configuration 
phi_v = fzero(@(phi) sincc(phi/pi)-(obj.L+stroke)/obj.lv,[0 pi]);
rv = obj.lv/(2*phi_v);

% find the position of the crest and the valley points in the deformed
% configuration
Rom = obj.Roe+rc*(1-cos(phi_c));
Rim = obj.Rie+rv*(1-cos(phi_v));

% find other geometric parameters of the deformed configuration
sol = obj.FindGeom(Rom,Rim,obj.alpha,obj.delta,obj.gamma,initial_guess);

Rcm = sol(1);
Rvm = sol(2);
beta_m = sol(3);

if beta_m>=pi || beta_m<=obj.alpha/2, StopFlag = true; end

% maximum latitude angle
vmax = atan((obj.L+stroke)/(2*obj.Roe));

rc0 = Rom-rc;

s1x_m = Rcm*sin(beta_m);
s1y_m = Rcm*(cos(beta_m)-1)+Rom;

os1_m = sqrt(s1x_m^2+s1y_m^2);
os2_m = Rvm*sin(beta_m-obj.alpha/2)/sin(obj.gamma);

[x0s1,~,rs1] = obj.CircleFrom3Points([obj.os1_e,-(obj.L+stroke)/2],[os1_m,0],[obj.os1_e,(obj.L+stroke)/2]);
[x0s2,~,rs2] = obj.CircleFrom3Points([obj.os2_e,-(obj.L+stroke)/2],[os2_m,0],[obj.os2_e,(obj.L+stroke)/2]);

phi_s1 = asin((obj.L+stroke)/(2*rs1));
phi_s2 = asin((obj.L+stroke)/(2*rs2));

if (phi_c>pi/2) || (phi_v>pi/2) || (phi_s1>pi/2) || (phi_s2>pi/2)
    StopFlag = true;
end

DeformedConfigParameters = [rc0,rc,phi_c,phi_v,rs1,x0s1,rs2,x0s2,vmax];

end