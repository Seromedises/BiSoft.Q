classdef Bellows < handle

    properties
        L 
        Rib 
        ns 
        alpha_s 
        Rs1,Rs2 

        d,D,h,lhs 
        AB,ltot 

        h1,Rt1,rcs1 
        h2,Rt2,rcs2 

        % parameters that change in deformed configuration
        DeformedConfigsParameters = struct('alpha_s_vec',[],'h_vec',[], ...
            'rcs1_vec',[],'Rs1_vec',[],'Rt1_vec',[],'h1_vec',[],'rcs2_vec',[], ...
            'Rs2_vec',[],'Rt2_vec',[],'h2_vec',[],'x_vec',[]);

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
        
        function obj = Bellows(L,Rib,ns,alpha_s,Rs1,Rs2)

            arguments
               L (1,1) double {mustBePositive}
               Rib (1,1) double {mustBePositive}
               ns (1,1) double {mustBeInteger,mustBePositive}
               alpha_s (1,1) double {mustBePositive,mustBeInRange(alpha_s,0,1.5708,"exclusive")} 
               Rs1 (1,1) double {mustBePositive} 
               Rs2 (1,1) double {mustBePositive} 
            end

            % Constructor
            obj.L = L;
            obj.Rib =Rib;
            obj.ns = ns;
            obj.alpha_s = alpha_s;
            obj.Rs1 = Rs1;
            obj.Rs2 = Rs2;

            % find the other parameters that define the geometry of the bellows
            obj.d = 2*(obj.Rib-obj.Rs1*(1-sin(obj.alpha_s))/sin(obj.alpha_s));
            obj.h = obj.L/(2*obj.ns);
            obj.lhs = obj.h/tan(obj.alpha_s);
            obj.D = obj.d+2*obj.lhs;
            
            obj.h1 = obj.Rs1*cos(obj.alpha_s);
            obj.h2 = obj.lhs*tan(obj.alpha_s)-obj.Rs2*cos(obj.alpha_s);
            
            obj.AB = obj.lhs/cos(obj.alpha_s);
            AA1 = obj.Rs1/tan(obj.alpha_s);
            A2B = obj.Rs2/tan(obj.alpha_s);
            
            obj.Rt1 = obj.d/2+AA1*cos(obj.alpha_s);
            obj.Rt2 = obj.d/2+(obj.AB-A2B)*cos(obj.alpha_s);
            
            obj.rcs1 = obj.Rib+obj.Rs1;
            obj.rcs2 = obj.d/2+obj.lhs-obj.Rs2/sin(obj.alpha_s);
            
            obj.ltot = obj.Rs1*(pi/2-obj.alpha_s)+(obj.h2-obj.h1)/sin(obj.alpha_s)+obj.Rs2*(pi/2-obj.alpha_s);
            
        end
        
        k_bell = AxialStiffness(obj,stroke)
        DeformedConfigParameters = DeformedConfiguration(obj,stroke)
        PlotNominalGeometry(obj)
        F_push = PushForce(obj,p)
        volume = Volume(obj,GeometricParams)

    end

    methods(Access=public,Hidden=true)

        [X,Y,Z] = BellowsSurface(obj,alpha_s,h,rcs1,Rs1,Rt1,h1,...
                  rcs2,Rs2,Rt2,stroke,section_view,angle1,angle2)

        [x_mirror,y_mirror] = MirrorPoints(obj,angle,xp,yp)

    end
    
end




