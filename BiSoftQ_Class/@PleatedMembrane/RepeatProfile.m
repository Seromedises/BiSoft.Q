function [x_profile,y_profile] = repeat_profile(~,x_curve,y_curve,n)

% INPUTS:
% - x_curve = column vector 
% - y_curve = column vector
% - n = number of times to repeat the curve

% OUTPUTS:
% - x_profile = column vector 
% - y_profile = column vector

rot_angle = linspace(0,2*pi,n);
x_profile = zeros(length(x_curve),length(rot_angle));
y_profile = zeros(length(y_curve),length(rot_angle));

for i=1:n
    Rot = [ cos(rot_angle(i)) sin(rot_angle(i));
           -sin(rot_angle(i)) cos(rot_angle(i))];
    profile = Rot*[x_curve y_curve]';
    x_profile(:,i) = profile(1,:)';
    y_profile(:,i) = profile(2,:)';
end

x_profile = reshape(x_profile,[],1);
y_profile = reshape(y_profile,[],1);

end