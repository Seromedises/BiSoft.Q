function [s1x,s1y,s2x,s2y] = InnerTangent(~,c1,r1,c2,r2)

% INPUTS:
% - c1 = vector of the coordinates of the center of circle 1 (dim 2x1)
% - r1 = radius of circle 1
% - c2 = vector of the coordinates of the center of circle 1 (dim 2x1)
% - r2 = radius of circle 2

% OUTPUTS:
% - s1x,s1y = x,y coordinates of the point of tangency on circle 1 
% - s2x,s2y = x,y coordinates of the point of tangency on circle 2

% N.B. the function outputs only the points of tangency of one of the two
% possible inner tangents (the one that produces a valid geometry)

xc1 = c1(1);
yc1 = c1(2);

xc2 = c2(1);
yc2 = c2(2);

hyp = sqrt((xc1-xc2)^2+(yc1-yc2)^2);
short = r1+r2;

angle = atan2(yc2-yc1,xc2-xc1)-asin(short/hyp)+pi/2;

s1x = xc1+r1*cos(angle);
s1y = yc1+r1*sin(angle);

s2x = xc2+r2*cos(angle+pi);
s2y = yc2+r2*sin(angle+pi);

if s2y > s1y
    % consider the other inner tangent. Otherwise the geometry will be
    % invalid
    angle = atan2(yc2-yc1,xc2-xc1)+asin(short/hyp)-pi/2;

    s1x = xc1+r1*cos(angle);
    s1y = yc1+r1*sin(angle);

    s2x = xc2+r2*cos(angle+pi);
    s2y = yc2+r2*sin(angle+pi);
end

end

