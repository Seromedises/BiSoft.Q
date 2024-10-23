classdef Bellows < handle

    properties
        L 
        Rib 
        nb 
        alpha_b 
        Rvb,Rcb 

        d,D,h,lhs 
        AB,ltot 

        h1,Rt1,rcvb 
        h2,Rt2,rccb 

        % parameters that change in deformed configuration
        DeformedConfigsParameters = struct('alpha_b_vec',[],'h_vec',[], ...
            'rcvb_vec',[],'Rvb_vec',[],'Rt1_vec',[],'h1_vec',[],'rccb_vec',[], ...
            'Rcb_vec',[],'Rt2_vec',[],'h2_vec',[],'x_vec',[]);

        % section view properties
        section_view = 0; % flag to enable section view when plotting the 3d geometry
        angle1_section_view = 0; % angle in radians
        angle2_section_view = 0; % angle in radians
        
        t = 1.2; % thickness of the membrane (mm)

        % material properties
        E = 10; % Young modulus (MPa)
        nu = 0.48; % Poisson coefficient

    end

    methods(Access=public)

        % TO DO: add method to validate geometry of the bellows
        
        function obj = Bellows(L,Rib,nb,alpha_b,Rvb,Rcb)

            arguments
               L (1,1) double {mustBePositive}
               Rib (1,1) double {mustBePositive}
               nb (1,1) double {mustBeInteger,mustBePositive}
               alpha_b (1,1) double {mustBePositive,mustBeInRange(alpha_b,0,1.5708,"exclusive")} 
               Rvb (1,1) double {mustBePositive} 
               Rcb (1,1) double {mustBePositive} 
            end

            % Constructor
            obj.L = L;
            obj.Rib =Rib;
            obj.nb = nb;
            obj.alpha_b = alpha_b;
            obj.Rvb = Rvb;
            obj.Rcb = Rcb;

            % find the other parameters that define the geometry of the bellows
            obj.d = 2*(obj.Rib-obj.Rvb*(1-sin(obj.alpha_b))/sin(obj.alpha_b));
            obj.h = obj.L/(2*obj.nb);
            obj.lhs = obj.h/tan(obj.alpha_b);
            obj.D = obj.d+2*obj.lhs;
            
            obj.h1 = obj.Rvb*cos(obj.alpha_b);
            obj.h2 = obj.lhs*tan(obj.alpha_b)-obj.Rcb*cos(obj.alpha_b);
            
            obj.AB = obj.lhs/cos(obj.alpha_b);
            AA1 = obj.Rvb/tan(obj.alpha_b);
            A2B = obj.Rcb/tan(obj.alpha_b);
            
            obj.Rt1 = obj.d/2+AA1*cos(obj.alpha_b);
            obj.Rt2 = obj.d/2+(obj.AB-A2B)*cos(obj.alpha_b);
            
            obj.rcvb = obj.Rib+obj.Rvb;
            obj.rccb = obj.d/2+obj.lhs-obj.Rcb/sin(obj.alpha_b);
            
            obj.ltot = obj.Rvb*(pi/2-obj.alpha_b)+(obj.h2-obj.h1)/sin(obj.alpha_b)+obj.Rcb*(pi/2-obj.alpha_b);
            
        end
        
        k_bell = AxialStiffness(obj,stroke)
        DeformedConfigParameters = DeformedConfiguration(obj,stroke)
        PlotNominalGeometry(obj)
        F_push = PushForce(obj,p)
        volume = Volume(obj,GeometricParams)

        TR = Surf2Solid(onj,offsetIn,offsetOut)

    end

    methods(Access=public,Hidden=true)

        [X,Y,Z,Nx,Ny,Nz] = BellowsSurface(obj,alpha_b,h,rcvb,Rvb,Rt1,h1,...
                  rccb,Rcb,Rt2,stroke,section_view,angle1,angle2)

        [x_mirror,y_mirror] = MirrorPoints(obj,angle,xp,yp)

    end
    
end




