classdef PleatedMembrane < handle

    properties
        % adimensional parameters
        a1,a2,a3,a4,a5,a6,a7  

        % number of pleats <-> angular distance between pleats
        Np  
        alphap 
        
        % longitudinal section parameters
        L 
        lcp,phi_cp,rc,rc0 
        lvp,phi_vp,rv 
        vmax 
        phi_s1,rs1,x0s1 
        phi_s2,rs2,x0s2 
        
        % cross section parameters
        Rib 
        Roep,Riep,Rcep,Rvep,beta_e,os1_e,os2_e 
        Romp,Rimp,Rcmp,Rvmp,beta_m 
        delta,gamma 
        area_ant 
 
        % theoretic strokes
        xtot_th 
        xmax_th 
        xmin_th 

        % parameter that change in deformed configuration
        DeformedConfigsParameters = struct('rc0_vec',[],'rc_vec',[],'phi_cp_vec',[], ...
            'phi_vp_vec',[],'rs1_vec',[],'x0s1_vec',[],'rs2_vec',[],...
            'x0s2_vec',[],'vmax',[],'x_vec',[]);
        
        ParametricStudyBool
        ValidGeometry

    end

    methods(Access=public)
        
        function obj = PleatedMembrane(a2,a3,a4,a5,a6,N,L,Riep,Rcep,Rvep,Rib,ParametricStudyBool)

            arguments
               a2 (1,1) double {mustBeInRange(a2,0.6366,1,"exclude-upper")} 
               a3 (1,1) double {mustBeGreaterThan(a3,0)}
               a4 (1,1) double {mustBeInRange(a4,0,1,"exclusive")}
               a5 (1,1) double {mustBeInRange(a5,0,1,"exclusive")}
               a6 (1,1) double {mustBeInRange(a6,a2,1)}
               N (1,1) double {mustBeInteger,mustBePositive}
               L (1,1) double {mustBePositive}
               Riep (1,1) double {mustBePositive}
               Rcep (1,1) double {mustBeNonzero}
               Rvep (1,1) double {mustBeNonzero}
               Rib (1,1) double {mustBePositive}
               ParametricStudyBool (1,1) logical = false
            end
            
            lastwarn("")
            warning('off')
            GeomDefWarnMsg = strings;

            % Constructor
            obj.a2 = a2;
            obj.a3 = a3;
            obj.a4 = a4;
            obj.a5 = a5;
            obj.a6 = a6;
            obj.Np = N;
            obj.L = L;
            obj.Riep = Riep;
            obj.Rcep = Rcep;
            obj.Rvep = Rvep;
            obj.Rib = Rib;
            obj.ParametricStudyBool = ParametricStudyBool;

            % find the other geometric parameters of the pleated membrane
            obj.lcp = obj.L/obj.a2;
            obj.lvp = obj.lcp*obj.a6;

            obj.a1 = obj.lvp/(2*obj.Riep);
            
            % compute theoric strokes
            obj.xmax_th = obj.lvp-obj.L; % maximum elongation stroke
            obj.xmin_th = obj.lcp*2/pi-obj.L; % maximum contraction stroke

            obj.xtot_th = obj.xmax_th-obj.xmin_th;

            % angle between two successive crests
            obj.alphap = 2*pi/obj.Np;

            obj.phi_vp = fzero(@(phi) sincc(phi/pi)-obj.L/obj.lvp,[0 pi]);
            obj.rv = obj.lvp/(2*obj.phi_vp);
            obj.Rimp = obj.Riep+obj.rv*(1-cos(obj.phi_vp));
            
            obj.Romp = obj.Rimp/obj.a5;
            obj.phi_cp = fzero(@(phi) sincc(phi/pi)-obj.L/obj.lcp,[0 pi]);
            obj.rc = obj.lcp/(2*obj.phi_cp);
            obj.rc0 = obj.Romp-obj.rc;
            obj.Roep = obj.Romp-obj.rc*(1-cos(obj.phi_cp));
             
            phi_c_xmin = fzero(@(phi) sincc(phi/pi)-(obj.L+obj.xmin_th)/obj.lcp,[0 pi]);
            rc_xmin = obj.lcp/(2*phi_c_xmin);
            Rom_xmin = obj.Roep+rc_xmin*(1-cos(phi_c_xmin));
            obj.a7 = Rom_xmin/obj.Rib;
            
            if isnan(obj.Rcep) || isnan(obj.Rvep)
                Rcep_max = fzero(@(Rcep_max) 1+a3-sqrt(((obj.Riep/Rcep_max+a3)*sin(obj.alphap/2))^2+(obj.Roep/Rcep_max-1-(obj.Riep/Rcep_max+obj.a3)*cos(obj.alphap/2))^2),1e-6);
                obj.Rcep = obj.a4*Rcep_max;
                obj.Rvep = obj.a3*obj.Rcep;
            end
            
            R1e = obj.Roep-obj.Rcep;
            R2e = obj.Riep+obj.Rvep;
            
            % centre of the crest fillet in the end plate cross section
            xc1_e = 0;
            yc1_e = R1e;
            
            % centre of the valley fillet in the end plate cross section
            xc2_e = R2e*sin(obj.alphap/2);
            yc2_e = R2e*cos(obj.alphap/2);
            
            dist_c1c2 = sqrt((xc1_e-xc2_e)^2+(yc1_e-yc2_e)^2);
            
            if obj.Rcep+obj.Rvep>=dist_c1c2
                % geometry not valid: the crest fillet and the
                % valley fillet are secant (>) or tangent (=) 
                warning("Geometry not valid (crest and valley fillet are secant or tangent)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end
            
            % find the inner tangent line and the points of tangency
            [s1x_e,s1y_e,s2x_e,s2y_e] = InnerTangent(obj,[xc1_e;yc1_e],obj.Rcep,[xc2_e;yc2_e],obj.Rvep);
            
            if (s2y_e>s1y_e) || (obj.Riep*sin(obj.alphap/2)<s2x_e) || (obj.Riep*cos(obj.alphap/2)>s2y_e) || (s1x_e<0) || (s1y_e<0) || (s2x_e<0) || (s2y_e<0) 
                % geometry not valid: the definition of crests and
                % valleys is not respected (first condition), or s2 is beyond the mirror line
                % passing through the valley (second and third
                % conditions),or the points of tangency are not in
                % the first quadrant (other conditions
                warning("Geometry not valid (check points of tangency)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end
            
            if ~isreal(s1x_e) || ~isreal(s1y_e) || ~isreal(s2x_e) || ~isreal(s2y_e) 
                % geometry not valid: points of tangency must be
                % real numbers
                warning("Geometry not valid (imaginary points of tangency)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end
            
            obj.os1_e = sqrt(s1x_e^2+s1y_e^2);
            obj.os2_e = sqrt(s2x_e^2+s2y_e^2);
            
            % beta angle on the end plate cross section
            obj.beta_e = pi-atan2(s1y_e-s2y_e,s1x_e-s2x_e);
            
            if obj.beta_e>=pi
                % geometry not valid: beta_e must be < pi 
                warning("Geometry not valid (beta_e>=pi)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end
                                 
            % maximum latitude angle in the nominal configuration
            obj.vmax = atan(obj.L/(2*obj.Roep));
            
            % define some useful adimensional parameters
            Ae = obj.Roep/obj.Rcep;
            Be = obj.Riep/obj.Rvep; 
            
            obj.delta = atan(sin(obj.beta_e)/(Ae-1+cos(obj.beta_e))); % the same for all cross sections
            obj.gamma = atan(sin(obj.beta_e-obj.alphap/2)/(Be+1-cos(obj.beta_e-obj.alphap/2))); % the same for all cross section
            
            if (obj.delta<0) || (obj.gamma<0) || (obj.vmax<0)
                % geometry not valid: these angles must be >=0
                warning("Geometry not valid (delta<0 or gamma<0 or vmax<0)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end
            
            % antagonistic area
            s1s2_e = sqrt((s1x_e-s2x_e)^2+(s1y_e-s2y_e)^2); % s1-s2 segment on the end plate 
            
            area1 = obj.Rcep^2*obj.beta_e/2;
            area2 = obj.Rcep*sin(obj.beta_e)*R1e/2;
            area3 = obj.os1_e*s1s2_e*sin(pi/2-(obj.beta_e-obj.delta))/2;
            area4 = R2e*obj.Rvep*sin(obj.beta_e-obj.alphap/2)/2-obj.Rvep^2*(obj.beta_e-obj.alphap/2)/2;
            
            if area1<0 || area2<0 || area3<0 || area4<0
                warning("Geometry not valid (area1, area2, area3 or area4 < 0)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end

            obj.area_ant = 2*N*(area1+area2+area3+area4);
            
            % find unknown geometric parameters of the middle cross section. We must
            % solve a non-linear system of three equations to find Rcmp, Rvmp, beta_m
            
            % N.B. the non-linear system of equations has several solutions. We must
            % find the solution for which the following condition is met: 
            % - delta and gamma on the middle cross section are equal to the ones on
            % the end plate cross section
            
            flag_geom = 0; 
            delta_m = 0; % initialize value of delta angle on the middle cross section
            gamma_m = 0; % initialize value of gamma angle on the middle cross section
            beta_guess_0 = obj.beta_e;
            beta_guess_1 = obj.beta_e;
            initial_guess = [obj.Romp/Ae; obj.Rimp/Be; beta_guess_0];
            phi_lb = 0; % initialize value of phi angle for u=0
            phi_ub = 0; % initialize value of phi angle for u=alpha/2
            
            while (abs(delta_m-obj.delta)>1e-6) || (abs(gamma_m-obj.gamma)>1e-6) || (abs(obj.phi_cp-phi_lb)>1e-6) || (abs(obj.phi_vp-phi_ub)>1e-6)
                
                sol = obj.FindGeom(obj.Romp,obj.Rimp,obj.alphap,obj.delta,obj.gamma,initial_guess);
            
                obj.Rcmp = sol(1);
                obj.Rvmp = sol(2);
                obj.beta_m = sol(3);
            
                if obj.beta_m<=obj.alphap/2 || obj.beta_m>=pi 
                    % geometry not valid
                    warning("Geometry not valid (beta_m<=alpha/2 or beta_m>=pi)")
                    if strcmp(GeomDefWarnMsg,"")
                        GeomDefWarnMsg = lastwarn;
                    end
                    break
                end
                
                % evaluate delta and gamma on the middle cross section
                delta_m = atan(sin(obj.beta_m)/(obj.Romp/obj.Rcmp-1+cos(obj.beta_m)));
                gamma_m = atan(sin(obj.beta_m-obj.alphap/2)/(obj.Rimp/obj.Rvmp+1-cos(obj.beta_m-obj.alphap/2)));
                   
                % evaluate phi for u=0 and u=alpha/2
                R1m = obj.Romp-obj.Rcmp;
                % R2m = obj.Rim+obj.Rvm;
            
                os1_m = sqrt(R1m^2+obj.Rcmp^2-2*R1m*obj.Rcmp*cos(pi-obj.beta_m)); % O-s1 segment on the middle cross section
                os2_m = obj.Rvmp*sin(obj.beta_m-obj.alphap/2)/sin(obj.gamma); % O-s2 segment on the middle cross section
                    
                [obj.x0s1,~,obj.rs1] = obj.CircleFrom3Points([obj.os1_e,-L/2],[os1_m,0],[obj.os1_e,L/2]);
                [obj.x0s2,~,obj.rs2] = obj.CircleFrom3Points([obj.os2_e,-L/2],[os2_m,0],[obj.os2_e,L/2]);
                
                [phi_func1,~,phi_func3] = MyFuncPhi(obj,obj.rc,obj.rc0,obj.rs1,obj.x0s1,obj.rs2,obj.x0s2,0,obj.vmax);
                
                phi_lb = phi_func1(0);
                phi_ub = phi_func3(obj.alphap/2);
            
                % update initial guess of beta_m 
                if flag_geom==0
                    beta_guess_0 = 0.99*beta_guess_0;
                    initial_guess = [obj.Romp/Ae; obj.Rimp/Be; beta_guess_0];
                elseif flag_geom==1
                    beta_guess_1 = 1.01*beta_guess_1;
                    initial_guess = [obj.Romp/Ae; obj.Rimp/Be; beta_guess_1];
                end
                
                if beta_guess_0<=obj.alphap/2
                    flag_geom = 1;
                end
                
                if flag_geom==1
                    if beta_guess_1>=pi
                        warning("It was not possible to find valid values for Rcm and Rvm")
                        if strcmp(GeomDefWarnMsg,"")
                            GeomDefWarnMsg = lastwarn;
                        end
                        break
                    end
                end
            
            end
            
            obj.phi_s1 = asin(obj.L/(2*obj.rs1));
            obj.phi_s2 = asin(obj.L/(2*obj.rs2));
            
            if (obj.phi_cp>pi/2) || (obj.phi_vp>pi/2) || (obj.phi_s1>pi/2) || (obj.phi_s2>pi/2)
                % geometry not valid: the angle is too big
                warning("Geometry not valid (phi_[]>pi/2)")
                if strcmp(GeomDefWarnMsg,"")
                    GeomDefWarnMsg = lastwarn;
                end
            end

            obj.ValidateGeometry(GeomDefWarnMsg);
            warning('on')

        end
            
        [DeformedConfigParameters,sol,StopFlag] = DeformedConfiguration(obj,...
                                                  stroke,initial_guess)

        F_pull = PullForce(obj,GeometricParams,p,stroke,method)
        F_pull_approx = PullForceApprox(obj,GeometricParams,p,stroke,method,approx_method)
        PlotNominalGeometry(obj)
        volume = Volume(obj,GeometricParams,stroke,ApproximatedBool)
        
        TR = Surf2Solid(obj,offsetIn,offsetOut)

    end
    
    methods(Access=public,Hidden=true)

        [X,Y,Z,X_curve,Y_curve,X_profile,Y_profile,Z_profile,Nx,Ny,Nz] = PleatedMembraneSurface(...
                                                        obj,rc0,rc,rs1,x0s1,rs2,x0s2,vmax)
        [x_mirror,y_mirror] = MirrorPoints(obj,angle,xp,yp)
    
    end
    
    methods(Access=protected,Hidden=true)
        
        [xc,yc,r] = CircleFrom3Points(obj,p1,p2,p3)
        F = FindGeom(obj,Ro,Ri,alpha,delta,gamma,x0)
        [s1x,s1y,s2x,s2y] = InnerTangent(obj,c1,r1,c2,r2)
        [phi_func1,phi_func2,phi_func3] = MyFuncPhi(obj,rc,rc0,rs1,x0s1,rs2,x0s2,...
                                                    stroke,vmax)

        [x_profile,y_profile] = RepeatProfile(obj,x_curve,y_curve,n)

        ValidateGeometry(obj,GeomDefWarnMsg)

    end

end

