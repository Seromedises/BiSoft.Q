function TR = ExportSTL(obj,offsetInPleat,offsetOutPleat,offsetInBell,offsetOutBell,...
    HeightCap,RadiusEndCap)

% Pleated Membrane
TR_Pleat = obj.PleatedMembrane.Surf2Solid(offsetInPleat,offsetOutPleat);
PleatedMembranePoints = TR_Pleat.Points;
PleatedMembraneFaces = TR_Pleat.ConnectivityList;

% Bellows
TR_Bell = obj.Bellows.Surf2Solid(offsetInBell,offsetOutBell);
BellowsPoints = TR_Bell.Points;
BellowsFaces = TR_Bell.ConnectivityList;

% Upper and Lower Cap
t_LowerCap = linspace(0,2*pi,40)';
x_LowerCap = RadiusEndCap.*cos(t_LowerCap);
y_LowerCap = RadiusEndCap.*sin(t_LowerCap);

z_LowerCap = ones(numel(x_LowerCap),1);
x_LowerCap = repmat(x_LowerCap,2,1);
y_LowerCap = repmat(y_LowerCap,2,1);
z_LowerCap = z_LowerCap*linspace(0,HeightCap,2);
z_LowerCap = z_LowerCap(:);
LowerCapShape = alphaShape(x_LowerCap,y_LowerCap,z_LowerCap);

LowerCapShape.Alpha = inf;
[elements,nodes] = boundaryFacets(LowerCapShape);
TR_LowerCap = triangulation(elements,nodes);

UpperCapPoints = TR_LowerCap.Points;
UpperCapPoints(:,3) = UpperCapPoints(:,3)+obj.PleatedMembrane.L/2;
UpperCapFaces = TR_LowerCap.ConnectivityList;

LowerCapPoints = TR_LowerCap.Points;
LowerCapPoints(:,3) = LowerCapPoints(:,3)-obj.PleatedMembrane.L/2-HeightCap;
LowerCapFaces = TR_LowerCap.ConnectivityList;

% Assemble BiSoft
allVertices = [PleatedMembranePoints; BellowsPoints; UpperCapPoints; LowerCapPoints];
allFaces = [PleatedMembraneFaces; BellowsFaces+size(PleatedMembranePoints,1); ...
    UpperCapFaces+size(PleatedMembranePoints,1)+size(BellowsPoints,1); ...
    LowerCapFaces+size(PleatedMembranePoints,1)+size(BellowsPoints,1)+size(UpperCapPoints,1)];

TR = triangulation(allFaces,allVertices);

figure()
trisurf(TR)
axis equal

stlwrite(TR,'monolithic_bisoft.stl')

end