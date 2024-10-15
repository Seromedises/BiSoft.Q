function [xc,yc,r] = CircleFrom3Points(~,p1,p2,p3)

% INPUTS:
% - p1 = vector of the (x,y) coordinates of the first point
% - p2 = vector of the (x,y) coordinates of the second point
% - p3 = vector of the (x,y) coordinates of the third point

% OUTPUTS:
% - xc,yc = (x,y) coordinates of the center of the circle
% - r = radius of the circle
    
z1 = p1(1)+1j*p1(2);
z2 = p2(1)+1j*p2(2);
z3 = p3(1)+1j*p3(2);

if (z1==z2) || (z2==z3) || (z3==z1)
    error("Duplicated points")
end
    
w = (z3-z1)/(z2-z1);

if abs(imag(w)) <= 0
    error("Points are collinear")
end
    
c = (z2-z1)*(w-abs(w)^2)/(2j*imag(w))+z1;
xc = real(c);
yc = imag(c);
r = abs(z1-c);

end
