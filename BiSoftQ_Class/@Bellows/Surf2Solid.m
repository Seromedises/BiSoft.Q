function TR = Surf2Solid(obj,offsetIn,offsetOut)

[X_profile,Y_profile,Z_profile,Nx,Ny,Nz] = ...
    obj.BellowsSurface(obj.alpha_b,obj.h,obj.rcvb,obj.Rvb,...
                       obj.Rt1,obj.h1,obj.rccb,obj.Rcb,obj.Rt2,0);

X_profile = X_profile';
Y_profile = Y_profile';
Z_profile = Z_profile';

Nx = Nx';
Ny = Ny';
Nz = Nz';

% unit normal vector
nx = Nx./sqrt(Nx.^2+Ny.^2+Nz.^2);
ny = Ny./sqrt(Nx.^2+Ny.^2+Nz.^2);
nz = Nz./sqrt(Nx.^2+Ny.^2+Nz.^2);

% offset surface - IN
Xt_in = X_profile-offsetIn*nx;
Yt_in = Y_profile-offsetIn*ny;
Zt_in = Z_profile-offsetIn*nz;

% offset surface - OUT
Xt_out = X_profile+offsetOut*nx;
Yt_out = Y_profile+offsetOut*ny;
Zt_out = Z_profile+offsetOut*nz;

% trim offset surface
[Zt_out(:,1),Zt_in(:,1)] = deal(Z_profile(1,1));
[Zt_out(:,end),Zt_in(:,end)] = deal(Z_profile(end,end));

[Ft_in,Vt_in] = mesh2tri(Xt_in,Yt_in,Zt_in,'b');
[Ft_out,Vt_out] = mesh2tri(Xt_out,Yt_out,Zt_out,'b');

Ft_in = fliplr(Ft_in);

V_wall_down = [[Xt_in(:,1); Xt_out(:,1)],[Yt_in(:,1); Yt_out(:,1)],[Zt_in(:,1); Zt_out(:,1)]];
V_wall_up = [[Xt_in(:,end); Xt_out(:,end)],[Yt_in(:,end); Yt_out(:,end)],[Zt_in(:,end); Zt_out(:,end)]]; 

nwv_down = length(V_wall_down)/2; % number of wall vertices
F_wall_down = zeros(2*(nwv_down-1),3);
for k = 1:nwv_down-1
    F_wall_down(k           ,:) = [k+1       , k      ,      k+nwv_down];
    F_wall_down(k+nwv_down-1,:) = [k+nwv_down, k+1+nwv_down, k+1];
end
% F_wall_down = fliplr(F_wall_down);

nwv_up = length(V_wall_up)/2; % number of wall vertices
F_wall_up = zeros(2*(nwv_up-1),3);
for k = 1:nwv_up-1
    F_wall_up(k         ,:) = [k+1     , k         , k+nwv_up];
    F_wall_up(k+nwv_up-1,:) = [k+nwv_up, k+1+nwv_up, k+1];
end

allVertices = [Vt_in; V_wall_down; Vt_out; V_wall_up];
allFaces = [Ft_in; F_wall_down+size(Vt_in,1); Ft_out+size(Vt_in,1)+size(V_wall_down,1); F_wall_up+size(Vt_in,1)+size(V_wall_down,1)+size(Vt_out,1)];

TR = triangulation(allFaces,allVertices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,V]=mesh2tri(X,Y,Z,tri_type)
% Available from http://www.mathworks.com/matlabcentral/fileexchange/28327
[J,I]=meshgrid(1:1:size(X,2)-1,1:1:size(X,1)-1);
switch tri_type
    case 'f' % Forward slash
        TRI_I=[I(:),I(:)+1,I(:)+1;  I(:),I(:),I(:)+1];
        TRI_J=[J(:),J(:)+1,J(:);   J(:),J(:)+1,J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'b' % Back slash
        TRI_I=[I(:),I(:)+1,I(:);  I(:)+1,I(:)+1,I(:)];
        TRI_J=[J(:)+1,J(:),J(:);   J(:)+1,J(:),J(:)+1];
        F = sub2ind(size(X),TRI_I,TRI_J);
    case 'x' % Cross
        TRI_I=[I(:)+1,I(:);  I(:)+1,I(:)+1;  I(:),I(:)+1;    I(:),I(:)];
        TRI_J=[J(:),J(:);    J(:)+1,J(:);    J(:)+1,J(:)+1;  J(:),J(:)+1];
        IND=((numel(X)+1):numel(X)+prod(size(X)-1))';
        F = sub2ind(size(X),TRI_I,TRI_J);
        F(:,3)=repmat(IND,[4,1]);
        Fe_I=[I(:),I(:)+1,I(:)+1,I(:)]; Fe_J=[J(:),J(:),J(:)+1,J(:)+1];
        Fe = sub2ind(size(X),Fe_I,Fe_J);
        Xe=mean(X(Fe),2); Ye=mean(Y(Fe),2);  Ze=mean(Z(Fe),2);
        X=[X(:);Xe(:)]; Y=[Y(:);Ye(:)]; Z=[Z(:);Ze(:)];
end
V=[X(:),Y(:),Z(:)];

end

end

