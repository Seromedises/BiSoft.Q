% This script shows the basic usage of the BiSoftQ class to design a
% BiSoft.Q actuator. A sample set of geometric parameters is provided for
% both type "A" and "B" actuators

% Matlab 2023b version required
% Parallel Computing Toolbox required
%

clear
close all
clc

addpath(genpath("BiSoftQ_Class")); 

tic 

run("SetPlotStyleProperties.m")

actuator_type = "A";

p_gauge_MPa = 0.08; % gauge pressure (MPa)

% geometric parameters of the pleated membrane
a2 = 0.85; % a2 = L/l_{c,p} (theoric values: 2/pi<=a2<l_{v,p}/l_{c,p})
a3 = 0.70; % a3 = R_{ve,p}/R_{ce,p}
a4 = 0.50; % a4 = R{ce,p}/R_{ce,p}max (0<a4<1)
a5 = 0.65; % a5 = R_{im,p}/R_{om,p} (<1)
a6 = 0.90; % a6 = l_{v,p}/l_{c,p} (<=1)

Np = 8;
L = 90;
Riep = 17;
Rcep = NaN;
Rvep = NaN;

% geometric parameters of the bellows
Rib = 11.25; % inner radius of the bellows 
nb = 5; % number of convolutions of the bellows
alpha_b = deg2rad(45); % half of the convolution angle
Rvb = 0.5; % fillet radius of the concave part
Rcb = 2; % fillet radius of the convex part

% membrane = PleatedMembrane(a2,a3,a4,a5,a6,Np,L,Riep,Rcep,Rvep,Rib);
% membrane.PlotNominalGeometry();
% bellows = Bellows(L,Rib,nb,alpha_b,Rvb,Rcb);
% bellows.PlotNominalGeometry();

NPoints = 30;

% Create an instance of BiSoftQ class
actuator = BiSoftQ(PleatedMembrane(a2,a3,a4,a5,a6,Np,L,Riep,Rcep,Rvep,Rib), ...
                   Bellows(L,Rib,nb,alpha_b,Rvb,Rcb), ...
                   actuator_type,p_gauge_MPa);

% actuator.Bellows.section_view = 1;
% actuator.Bellows.angle1_section_view = deg2rad(40);
% actuator.Bellows.angle2_section_view = deg2rad(220);

% Plot nominal geometry (3d surface, cross sections, longitudinal cross
% sections), i.e. geometry a x=0
actuator.PlotNominalGeometry();

% Show summary of design parameters. It is possible to save them in a .txt
% file
actuator.SummaryDesignParameters("SaveSummary",false,...
    "SummaryName","BiSoftQ_Design_Parameters_Summary.txt");

% Evaluate force characteristic of the actuator. Three different methods of
% evaluating the pulling force can be used. There are also two approximated
% methods. It is also possible to run this method as 'mexfile' (compiled,
% faster), or 'mfile'

% material properties (TPU Filaflex 82A)
actuator.Bellows.E = 22; % MPa, from datasheet
actuator.Bellows.nu = 0.39;
actuator.Bellows.t = 1; % mm

actuator.ForceCharacteristic("NPoints",NPoints,"Method","method_2","ApproxMethod","",...
        "FileType","mexfile","ElasticForceBool",true,"NumWorkers",0);

% Plot the force characteristic of the actuator
actuator.PlotForceCharacteristic();

S_pull_vec = (actuator.F_pull_vec-actuator.F_el_vec)/p_gauge_MPa; % mm^2
S_push_vec = (actuator.F_push_vec-actuator.F_el_vec)/p_gauge_MPa; % mm^2

% Evaluate and plot the volume of the chambers of the actuator
%Change NumWorkers (ex. inf) if you'd like to use the parallel computing
%functionality
actuator.VolumeChambers("NPoints",NPoints,"ApproxVolumeBool",false,"NumWorkers",0);
actuator.PlotVolumeChambers();

actuator.Efficiency();

% Create an animation of the actuator. It is possibile to save it as a
% video
actuator.AnimateActuator("NPoints",NPoints,"SaveVideo",false,...
        "VideoName","BiSoftQ_Animation","SeparateVideos",false,"NumWorkers",0)

actuator.ExportSTL(0.5,0.5,0.75,0.25,5,24.2);

toc
