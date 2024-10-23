function [] = ForceCharacteristic(obj,options)

% Method to compute the force characteristic of BiSoftQ ({F_pull,F_push} vs stroke)
%
% Several options can be used:
%
% - "Method": specify which method you want to use to compute
%   pulling force. The following values are possible:
%   * "method_1": pressure integral on the surface of the pleated membrane.
%   We assume that the resulting radial force is applied on the valley
%   fibers
%   * "method_2" (default): pressure integral on the surface of the pleated
%   membrane. We assume that the resulting radial force is applied on the
%   crest fibers
%   * "method_3": distributed force, integrated along the path described by
%   the pleated membrane on the terminal cross section
%
% - "ApproxMethod": specify which approximated method you want to use to
%   compute pulling force. The following values are possible:
%   * "" (default): do not approximate radial force
%   * "approx_1": radial force is computed as the pressure force acting on
%   a rectangle
%   * "approx_2": radial force is computed as the pressure force acting on
%   two trapezoids
%
% - "FileType": run this method either as a .m file, or as a .mex file.
%   Accepted values:
%   * "mfile" (default)
%   * "mexfile"
% 
% - "ElasticForceBool" (default=false): include or not an approximation of
%   the elastic force
% 
% - "NumWorkers" (default=0): number of workers to use in parfor loops. Set
%   to 0 to disable parfor loops (they will run as normal loops)

arguments
    obj BiSoftQ 
    options.NPoints (1,1) double = 100
    options.Method (1,1) string {mustBeMember(options.Method,{'method_1','method_2','method_3'})} = "method_2"
    options.ApproxMethod (1,1) string {mustBeMember(options.ApproxMethod,{'','approx_1','approx_2'})} = ""
    options.FileType (1,1) string {mustBeMember(options.FileType,{'mfile','mexfile'})} = "mfile"
    options.ElasticForceBool (1,1) logical = false
    options.NumWorkers (1,1) double = 0
end

NPoints = options.NPoints;
Method = options.Method;
ApproxMethod = options.ApproxMethod;
FileType = options.FileType;
ElasticForceBool = options.ElasticForceBool;
NumWorkers = options.NumWorkers;

DeformedConfigsParametersPlM = NaN(NPoints+1,2,9);
% order of parameters along third dimension of 'DeformedConfigsParameters':
% rc0, rc, phi_c, phi_v, rs1, x0s1, rs2, x0s2, vmax

DeformedConfigsParametersBel = NaN(NPoints+1,2,10);
% order of parameters along third dimension of 'DeformedConfigsParametersBel':
% alpha_s, h, rcs1, Rs1, Rt1, h1, rcs2, Rs2, Rt2, h2

if not(isempty(obj.F_pull_vec))  
    return
else
    x_vec = NaN(NPoints+1,2);
    x_vec(1,:) = 0;
    dx = 1/NPoints*[obj.PleatedMembrane.xmin_th; 0.95*obj.PleatedMembrane.xmax_th];
    obj.NPoints = NPoints;
end

NominalConfigParametersPlM = [obj.PleatedMembrane.rc0,obj.PleatedMembrane.rc,obj.PleatedMembrane.phi_cp,...
    obj.PleatedMembrane.phi_vp,obj.PleatedMembrane.rs1,obj.PleatedMembrane.x0s1,...
    obj.PleatedMembrane.rs2,obj.PleatedMembrane.x0s2,obj.PleatedMembrane.vmax];

NominalConfigParametersBel = [obj.Bellows.alpha_b,obj.Bellows.h,obj.Bellows.rcvb,...
    obj.Bellows.Rvb,obj.Bellows.Rt1,obj.Bellows.h1,obj.Bellows.rccb,obj.Bellows.Rcb,...
    obj.Bellows.Rt2,obj.Bellows.h2];

for jj=1:2
    DeformedConfigsParametersPlM(1,jj,:) = NominalConfigParametersPlM;
    DeformedConfigsParametersBel(1,jj,:) = NominalConfigParametersBel;
end

ConstantParametersPlM = [obj.PleatedMembrane.L,obj.PleatedMembrane.lcp,...
obj.PleatedMembrane.lvp,obj.PleatedMembrane.Roep,obj.PleatedMembrane.Riep,...
obj.PleatedMembrane.alphap,obj.PleatedMembrane.delta,obj.PleatedMembrane.gamma,...
obj.PleatedMembrane.os1_e,obj.PleatedMembrane.os2_e,obj.PleatedMembrane.area_ant];

F_pull_vec = NaN(NPoints+1,2);

if ~strcmp(ApproxMethod,"") 
    F_pull_vec(1,:) = obj.PleatedMembrane.PullForceApprox(NominalConfigParametersPlM,obj.p_gauge_MPa,0,Method,ApproxMethod);
elseif strcmp(FileType,"mexfile")
    F_pull_vec(1,:) = PullForce_mex(ConstantParametersPlM,NominalConfigParametersPlM,obj.p_gauge_MPa,0,Method);
else
    F_pull_vec(1,:) = obj.PleatedMembrane.PullForce(NominalConfigParametersPlM,obj.p_gauge_MPa,0,Method);
end

if ElasticForceBool
    k_bell_vec = NaN(NPoints+1,2);
    k_bell_vec(1,:) = obj.Bellows.AxialStiffness(0);
end

Rcmp = obj.PleatedMembrane.Rcmp;
Rvmp = obj.PleatedMembrane.Rvmp;
beta_m = obj.PleatedMembrane.beta_m;

iter_max = size(x_vec,1);

parfor (jj=1:2,NumWorkers) 

    initial_guess = [Rcmp,Rvmp,beta_m];
    
    for ii=2:iter_max
        
        x_vec(ii,jj) = dx(jj)*(ii-1);

        [DeformedConfigsParametersPlM(ii,jj,:),initial_guess,StopFlag] = ...
        obj.PleatedMembrane.DeformedConfiguration(x_vec(ii,jj),initial_guess);

        if StopFlag==true
            x_vec(ii,jj) = NaN;
            DeformedConfigsParametersPlM(ii,jj,:) = NaN;
            break
        end
        
        if  ~strcmp(ApproxMethod,"") 
            F_pull_vec(ii,jj) = obj.PleatedMembrane.PullForceApprox(squeeze(DeformedConfigsParametersPlM(ii,jj,:))',obj.p_gauge_MPa,x_vec(ii,jj),Method,ApproxMethod);
        elseif strcmp(FileType,"mexfile")
            F_pull_vec(ii,jj) = PullForce_mex(ConstantParametersPlM,squeeze(DeformedConfigsParametersPlM(ii,jj,:))',obj.p_gauge_MPa,x_vec(ii,jj),Method);
        else
            F_pull_vec(ii,jj) = obj.PleatedMembrane.PullForce(DeformedConfigsParametersPlM(ii,jj,:),obj.p_gauge_MPa,x_vec(ii,jj),Method);
        end
        
        if F_pull_vec(ii,jj)<0
            F_pull_vec(ii,jj) = NaN;
            x_vec(ii,jj) = NaN;
            DeformedConfigsParametersPlM(ii,jj,:);
            break
        end
        
        if ElasticForceBool
            k_bell_vec(ii,jj) = obj.Bellows.AxialStiffness(x_vec(ii,jj));
        end

        DeformedConfigsParametersBel(ii,jj,:) = obj.Bellows.DeformedConfiguration(x_vec(ii,jj));

    end

end

fns = fieldnames(obj.PleatedMembrane.DeformedConfigsParameters);
for idx=1:numel(fns)-1
    obj.PleatedMembrane.DeformedConfigsParameters.(fns{idx}) = DeformedConfigsParametersPlM(:,:,idx);
end
obj.PleatedMembrane.DeformedConfigsParameters.x_vec = x_vec;    
x_vec(:,1) = flip(x_vec(:,1));
x_vec = reshape(x_vec,[],1);
x_vec(NPoints+1) = [];

fns = fieldnames(obj.Bellows.DeformedConfigsParameters);
for idx=1:numel(fns)-1
   obj.Bellows.DeformedConfigsParameters.(fns{idx}) = DeformedConfigsParametersBel(:,:,idx);
end
obj.Bellows.DeformedConfigsParameters.x_vec = x_vec;

F_pull_vec(:,1) = flip(F_pull_vec(:,1));
F_pull_vec = reshape(F_pull_vec,[],1);
F_pull_vec(NPoints+1) = [];

if strcmp(obj.type,"A")
    F_pull_vec = F_pull_vec-obj.p_gauge_MPa*(obj.PleatedMembrane.area_ant-pi*obj.Bellows.Rib^2);
elseif strcmp(obj.type,"B")
    F_pull_vec = F_pull_vec-obj.p_gauge_MPa*obj.PleatedMembrane.area_ant;
end

if ElasticForceBool
    k_bell_vec(:,1) = flip(k_bell_vec(:,1));
    k_bell_vec = reshape(k_bell_vec,[],1);
    k_bell_vec(NPoints+1) = [];
    F_el_bell_vec = cumtrapz(x_vec,k_bell_vec);
    F_el_bell_vec = F_el_bell_vec-F_el_bell_vec(x_vec==0);
    F_el_membr_vec = zeros(size(F_el_bell_vec)); % TO DO: model stiffness of the pleated membrane
else
    [F_el_membr_vec,F_el_bell_vec] = deal(zeros(size(x_vec)));
end

F_el = F_el_membr_vec+F_el_bell_vec;

idx_ltz = find(F_pull_vec+F_el<=0,1,"last");
if ~isempty(idx_ltz)   
    [F_pull_vec(1:idx_ltz),x_vec(1:idx_ltz),F_el_membr_vec(1:idx_ltz),F_el_bell_vec(1:idx_ltz)] = deal(NaN);
end

idx_rtz = find(obj.Bellows.PushForce(obj.p_gauge_MPa)*ones(size(x_vec))+F_el>=0,1,"first");
if ~isempty(idx_rtz)   
    [F_pull_vec(idx_rtz:end),x_vec(idx_rtz:end),F_el_membr_vec(idx_rtz:end),F_el_bell_vec(idx_rtz:end)] = deal(NaN);
end

obj.F_el_vec = F_el_membr_vec+F_el_bell_vec;
obj.F_el_membr_vec = F_el_membr_vec; 
obj.F_el_bell_vec = F_el_bell_vec;

obj.x_vec = x_vec;
obj.F_pull_vec = F_pull_vec+obj.F_el_vec;
obj.F_push_vec = obj.Bellows.PushForce(obj.p_gauge_MPa)*ones(size(obj.x_vec))+obj.F_el_vec;

obj.xmin_eff = obj.x_vec(find(~isnan(obj.x_vec),1,'first'));
obj.xmax_eff = obj.x_vec(find(~isnan(obj.x_vec),1,'last'));
obj.xtot_eff = obj.xmax_eff-obj.xmin_eff;

end
