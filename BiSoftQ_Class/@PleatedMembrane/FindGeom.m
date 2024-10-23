% Notation:
% x(1) = Rc;
% x(2) = Rv;
% x(3) = beta;

% F(1) = 0 comes from the definition of delta angle
% F(2) = 0 comes from the definition of gamma angle
% F(3) = 0 comes from the formula for the distance of the point c2 from the
% line s1s2 (line s1s2 defined by point s1 and beta angle)

function F = FindGeom(~,Ro,Ri,alpha,delta,gamma,x0)

options = optimset('Display','off');
F = fsolve(@equations,x0,options);

    function F = equations(x)
        F(1) = Ro-x(1)*(1+sin(x(3))/tan(delta)-cos(x(3))); 
        F(2) = Ri-x(2)*(sin(x(3)-alpha/2)/tan(gamma)+cos(x(3)-alpha/2)-1); 
        F(3) = abs(x(1)/x(2)*(cos(x(3))*(1-Ro/x(1))-1)+(Ri/x(2)+1)*cos(x(3)-alpha/2))-1; 
    end

end