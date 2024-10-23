function  VolumeChambers(obj,options)

arguments
    obj BiSoftQ 
    options.NPoints (1,1) double = 100
    options.ApproxVolumeBool (1,1) logical = false
    options.NumWorkers (1,1) double = 0
end

NPoints = options.NPoints;
ApproxVolumeBool = options.ApproxVolumeBool;
NumWorkers = options.NumWorkers;

DeformedConfigsParametersPlM = NaN(NPoints+1,2,9);
% order of parameters along third dimension of 'DeformedConfigsParametersPlM':
% rc0, rc, phi_cp, phi_vp, rs1, x0s1, rs2, x0s2, vmax

DeformedConfigsParametersBel = NaN(NPoints+1,2,10);
% order of parameters along third dimension of 'DeformedConfigsParametersBel':
% alpha_s, h, rcs1, Rs1, Rt1, h1, rcs2, Rs2, Rt2, h2

if isempty(obj.F_pull_vec)
    % compute force characteristic
    obj.ForceCharacteristic("NPoints",NPoints,"Method","method_2","ApproxMethod","",...
        "FileType","mexfile","ElasticForceBool",false,"NumWorkers",0);
end

fns = fieldnames(obj.PleatedMembrane.DeformedConfigsParameters);
for idx=1:numel(fns)-1
   DeformedConfigsParametersPlM(:,:,idx) = obj.PleatedMembrane.DeformedConfigsParameters.(fns{idx});
end 
x_vec = obj.PleatedMembrane.DeformedConfigsParameters.x_vec;
idxs = zeros(2,1);
for ii=1:2
    if isempty(find(isnan(x_vec(:,ii)),1,'first'))
        idxs(ii) = length(x_vec(:,ii))+1;
    else
        idxs(ii) = find(isnan(x_vec(:,ii)),1,'first');
    end
end

fns = fieldnames(obj.Bellows.DeformedConfigsParameters);
for idx=1:numel(fns)-1
   DeformedConfigsParametersBel(:,:,idx) = obj.Bellows.DeformedConfigsParameters.(fns{idx});
end 

[volume_membr_vec,volume_bell_vec] = deal(NaN(NPoints+1,2));

NominalConfigParametersPlM = [obj.PleatedMembrane.rc0,obj.PleatedMembrane.rc,obj.PleatedMembrane.phi_cp,...
    obj.PleatedMembrane.phi_vp,obj.PleatedMembrane.rs1,obj.PleatedMembrane.x0s1,...
    obj.PleatedMembrane.rs2,obj.PleatedMembrane.x0s2,obj.PleatedMembrane.vmax];

volume_membr_vec(1,:) = obj.PleatedMembrane.Volume(NominalConfigParametersPlM,0,ApproxVolumeBool);

NominalConfigParametersBel = [obj.Bellows.alpha_b,obj.Bellows.h,obj.Bellows.rcvb,...
    obj.Bellows.Rvb,obj.Bellows.Rt1,obj.Bellows.h1,obj.Bellows.rccb,obj.Bellows.Rcb,...
    obj.Bellows.Rt2,obj.Bellows.h2];

volume_bell_vec(1,:) = obj.Bellows.Volume(NominalConfigParametersBel);

for jj=1:2
    DeformedConfigsParametersPlM(1,jj,:) = NominalConfigParametersPlM;
    DeformedConfigsParametersBel(1,jj,:) = NominalConfigParametersBel;
end

iter_max = size(x_vec,1);

parfor (jj=1:2,NumWorkers) 

    for ii=2:iter_max

        if ii==idxs(jj)
            break 
        else
            volume_membr_vec(ii,jj) = obj.PleatedMembrane.Volume(DeformedConfigsParametersPlM(ii,jj,:),x_vec(ii,jj),ApproxVolumeBool);
            volume_bell_vec(ii,jj) = obj.Bellows.Volume(DeformedConfigsParametersBel(ii,jj,:));
        end  

    end

end

volume_membr_vec(:,1) = flip(volume_membr_vec(:,1));
volume_membr_vec = reshape(volume_membr_vec,[],1);
volume_membr_vec(NPoints+1) = [];

volume_bell_vec(:,1) = flip(volume_bell_vec(:,1));
volume_bell_vec = reshape(volume_bell_vec,[],1);
volume_bell_vec(NPoints+1) = [];

if strcmp(obj.type,"A")

    obj.VolumeChamber1_vec = volume_membr_vec-volume_bell_vec;
    obj.VolumeChamber2_vec = volume_bell_vec;

elseif strcmp(obj.type,"B")

    obj.VolumeChamber1_vec = volume_membr_vec;
    obj.VolumeChamber2_vec = volume_bell_vec-volume_membr_vec;

end


end