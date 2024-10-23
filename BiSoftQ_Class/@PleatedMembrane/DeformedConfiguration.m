function [DeformedConfigParameters,sol,StopFlag] = DeformedConfiguration(obj,stroke,initial_guess)

StopFlag = false;

% find angle and radius of the crest arc in the deformed configuration 
phi_cp = fzero(@(phi) sincc(phi/pi)-(obj.L+stroke)/obj.lcp,[0 pi]);
rc = obj.lcp/(2*phi_cp);

% find angle and radius of the valley arc in the deformed configuration 
phi_vp = fzero(@(phi) sincc(phi/pi)-(obj.L+stroke)/obj.lvp,[0 pi]);
rv = obj.lvp/(2*phi_vp);

% find the position of the crest and the valley points in the deformed
% configuration
Romp = obj.Roep+rc*(1-cos(phi_cp));
Rimp = obj.Riep+rv*(1-cos(phi_vp));

% find other geometric parameters of the deformed configuration
sol = obj.FindGeom(Romp,Rimp,obj.alphap,obj.delta,obj.gamma,initial_guess);

Rcmp = sol(1);
Rvmp = sol(2);
beta_m = sol(3);

if beta_m>=pi || beta_m<=obj.alphap/2, StopFlag = true; end

% maximum latitude angle
vmax = atan((obj.L+stroke)/(2*obj.Roep));

rc0 = Romp-rc;

s1x_m = Rcmp*sin(beta_m);
s1y_m = Rcmp*(cos(beta_m)-1)+Romp;

os1_m = sqrt(s1x_m^2+s1y_m^2);
os2_m = Rvmp*sin(beta_m-obj.alphap/2)/sin(obj.gamma);

[x0s1,~,rs1] = obj.CircleFrom3Points([obj.os1_e,-(obj.L+stroke)/2],[os1_m,0],[obj.os1_e,(obj.L+stroke)/2]);
[x0s2,~,rs2] = obj.CircleFrom3Points([obj.os2_e,-(obj.L+stroke)/2],[os2_m,0],[obj.os2_e,(obj.L+stroke)/2]);

phi_s1 = asin((obj.L+stroke)/(2*rs1));
phi_s2 = asin((obj.L+stroke)/(2*rs2));

if (phi_cp>pi/2) || (phi_vp>pi/2) || (phi_s1>pi/2) || (phi_s2>pi/2)
    StopFlag = true;
end

DeformedConfigParameters = [rc0,rc,phi_cp,phi_vp,rs1,x0s1,rs2,x0s2,vmax];

end