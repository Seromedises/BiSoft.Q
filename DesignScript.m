% This script shows the basic usage of the BiSoftQ class to design a
% BiSoft.Q actuator. A sample set of geometric parameters is provided for
% both type "A" and "B" actuators

% Matlab 2023b version required

clear
close all
clc

addpath(genpath("BiSoftQ_Class")); % Remove it and add to path manually if you want to evaluate the computational cost
run("SetPlotStyleProperties.m")

tic 

% example of type "A" parameters
% ----------------------------------------------------------------------- %
% 
% actuator_type = "A";
% 
% p_gauge_MPa = 0.1; % gauge pressure (MPa)
% 
% % geometric parameters of the pleated membrane
% a2 = 0.85; % a2 = L/lc (theoric values: 2/pi<=a2<lv/lc)
% a3 = 0.70; % a3 = Rve/Rce
% a4 = 0.50; % a4 = Rce/Rce_max (0<a4<1)
% a5 = 0.70; % a5 = Rim/Rom (<1)
% a6 = 0.90; % a6 = lv/lc (<=1)
% 
% N = 8;
% L = 30;
% Rie = 15.48;
% Rce = NaN;
% Rve = NaN;
% 
% % geometric parameters of the bellows
% Rib = 11.10; % inner radius of the bellows
% ns = 3; % number of convolutions of the bellows
% alpha_s = deg2rad(45); % half of the convolution angle
% Rs1 = 1; % fillet radius of the concave part
% Rs2 = 3; % fillet radius of the convex part
% ----------------------------------------------------------------------- %


% example of type "B" parameters
% ----------------------------------------------------------------------- %

actuator_type = "B"; 

p_gauge_MPa = 0.1; % gauge pressure (MPa)

% geometric parameters of the pleated membrane
a2 = 0.85; % a2 = L/lc (theoric values: 2/pi<=a2<=lv/lc)
a3 = 0.70; % a3 = Rve/Rce
a4 = 0.98; % a4 = Rce/Rce_max (0<a4<1)
a5 = 0.658868116783836; % a5 = Rim/Rom (<1)
a6 = 0.90; % a6 = lv/lc (<=1)

N = 8;
L = 30;
Rie = 3.5;
Rce = NaN;
Rve = NaN;

% geometric parameters of the bellows
Rib = 19.414214; % inner radius of the bellows
ns = 3; % number of convolutions of the bellows
alpha_s = deg2rad(45); % half of the convolution angle
Rs1 = 1; % fillet radius of the concave part
Rs2 = 3; % fillet radius of the convex part
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% | The following part is the same for both type "A" and "B" actuators  | %
% ----------------------------------------------------------------------- %

% membrane = PleatedMembrane(a2,a3,a4,a5,a6,N,L,Rie,Rce,Rve,Rib);
% membrane.PlotNominalGeometry();
% bellows = Bellows(L,Rib,ns,alpha_s,Rs1,Rs2);
% bellows.PlotNominalGeometry();

NPoints = 100;

% Create an instance of BiSoftQ class
actuator = BiSoftQ(PleatedMembrane(a2,a3,a4,a5,a6,N,L,Rie,Rce,Rve,Rib), ...
                   Bellows(L,Rib,ns,alpha_s,Rs1,Rs2), ...
                   actuator_type,p_gauge_MPa);

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
%
% actuator.ForceCharacteristic("NPoints",NPoints,"Method","method_2","ApproxMethod","",...
%         "FileType","mexfile","ElasticForceBool",false,"NumWorkers",inf);

% If you want to include an approximation of the elastic force in the
% actuation forces, please run the following lines instead. TBN: at the
% moment there is no analytic model to evaluate the axial stiffness of the
% pleated membrane, thus the axial stiffness of this element is set equal
% to the axial stiffness of the bellows

% material properties (TPU 60A - TPU 82A)
actuator.Bellows.E = 10; % MPa
actuator.Bellows.nu = 0.48;
actuator.Bellows.t = 1.2; % mm

actuator.ForceCharacteristic("NPoints",NPoints,"Method","method_2","ApproxMethod","",...
        "FileType","mexfile","ElasticForceBool",true,"NumWorkers",inf);

% Plot the force characteristic of the actuator
actuator.PlotForceCharacteristic();

% S_pull_vec = actuator.F_pull_vec/p_gauge_MPa; % mm^2
% S_push_vec = actuator.F_push_vec/p_gauge_MPa; % mm^2

% Evaluate and plot the volume of the chambers of the actuator
actuator.VolumeChambers("NPoints",NPoints,"ApproxVolumeBool",true,"NumWorkers",inf);
actuator.PlotVolumeChambers();

% Create an animation of the actuator. It is possibile to save it as a
% video
% actuator.AnimateActuator("NPoints",NPoints,"SaveVideo",false,...
%         "VideoName","BiSoftQ_Animation","SeparateVideos",false)

toc



