function  VolumeChambers(obj)

arguments
    obj BiSoftQ
end

if isempty(obj.F_pull_vec)
    % compute force characteristic with an approx method
    obj.ForceCharacteristic("NPoints",NPoints,"method","method_2","approx_method","",...
        "filetype","mfile","ElasticForceBool",true,"NumWorkers",0);
end

% input work: \int_{in}^{out} Vdp (work required to compress the air)
% hypothesis: adiabatic compression, followed by isobaric cooling

Tamb = 293.15; % initial absolute temperature = final absolute temperature (K)
R = 287.1; % elastic constant of air (J/kg*K)
k = 1.4; % k = cp/cv
Pamb = 101325; % initial pressure in chamber 1 (at x=xmax), (Pa)

Ps = Pamb+(obj.p_gauge_MPa.*1e6); % final pressure in chamber 1 (at x=xmin), (Pa)

VolumeChamber1_vec = obj.VolumeChamber1_vec .* 1e-9;
VolumeChamber2_vec = obj.VolumeChamber2_vec .* 1e-9;

isnanindex = find(~isnan(obj.x_vec), 1);

if strcmp(obj.type,"A")

    deltaM_ch1 = (Ps .* VolumeChamber1_vec(isnanindex) - Pamb .* VolumeChamber1_vec(end)) ./ (R .* Tamb);
    deltaM_ch2 = (Ps .* VolumeChamber2_vec(end) - Pamb .* VolumeChamber2_vec(isnanindex)) ./ (R .* Tamb);
    deltaM = deltaM_ch1 + deltaM_ch2;
    
elseif strcmp(obj.type,"B")

    deltaM_ch1 = ((Ps .* VolumeChamber1_vec(isnanindex)) - (Pamb .* VolumeChamber1_vec(end))) ./ (R .* Tamb);
    deltaM_ch1_ch2 = ((Ps .* (VolumeChamber2_vec(end) + VolumeChamber1_vec(end))) - (Pamb .* (VolumeChamber2_vec(isnanindex) + VolumeChamber1_vec(isnanindex)))) ./ (R .* Tamb);
    deltaM = deltaM_ch1 + deltaM_ch1_ch2;
    
end

% input work
win = R .* Tamb .* (k / (k-1)) .* ((Ps/Pamb)^((k-1)/k) - 1); % J /kg Ideal work of compression for unit mass
Win = deltaM .* win; % J

% output work: \int_0^x F(x)dx
Wout_contr = cumtrapz(obj.x_vec(isnanindex:end)*1e-3,obj.F_pull_vec(isnanindex:end)); 
Wout_elo = abs(cumtrapz(obj.x_vec(isnanindex:end)*1e-3,obj.F_push_vec(isnanindex:end))); 
Wout = Wout_contr + Wout_elo;

obj.deltaM = deltaM; 
obj.Win = Win;
obj.Wout = Wout;

end