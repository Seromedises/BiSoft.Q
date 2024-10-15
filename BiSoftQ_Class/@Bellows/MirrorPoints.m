function [x_mirror,y_mirror] = MirrorPoints(~,angle,xp,yp)

% INPUTS:
% - angle = angle of the mirror line wrt x-axis, positive counterclockwise
% - xp = x-coordinates of the points to mirror
% - yp = y-coordinates of the points to mirror

% OUTPUTS:
% - x_mirror, y_mirror = x,y coordinates of the mirrored points

x_mirror = xp*cos(2*angle)-(-yp)*sin(2*angle);
y_mirror = xp*sin(2*angle)+(-yp)*cos(2*angle);

x_mirror = flip(x_mirror);
y_mirror = flip(y_mirror);

end
