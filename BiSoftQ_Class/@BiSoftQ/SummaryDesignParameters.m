function SummaryDesignParameters(obj,options)

arguments
    obj BiSoftQ 
    options.SaveSummary (1,1) logical = false
    options.SummaryName (1,:) string = "BiSoftQ_Design_Parameters_Summary.txt"
end

SaveSummary = options.SaveSummary;
SummaryName = options.SummaryName;

clc

if SaveSummary == 1, diary(SummaryName), end

fprintf('Summary of design parameters - pleated membrane \n')
fprintf('-----------------------------------------------\n')
fprintf('N     [-]: %i\n', obj.PleatedMembrane.N);
fprintf('xtot [mm]: %.6f\n', obj.PleatedMembrane.xtot_th);
fprintf('lc   [mm]: %.6f\n', obj.PleatedMembrane.lc);
fprintf('lv   [mm]: %.6f\n', obj.PleatedMembrane.lv);
fprintf('L    [mm]: %.6f\n', obj.PleatedMembrane.L);
fprintf('Roe  [mm]: %.6f\n', obj.PleatedMembrane.Roe);
fprintf('Rie  [mm]: %.6f\n', obj.PleatedMembrane.Rie);
fprintf('Rce  [mm]: %.6f\n', obj.PleatedMembrane.Rce);
fprintf('Rve  [mm]: %.6f\n', obj.PleatedMembrane.Rve);
fprintf('Rom  [mm]: %.6f\n', obj.PleatedMembrane.Rom);
fprintf('Rim  [mm]: %.6f\n', obj.PleatedMembrane.Rim);
fprintf('Rcm  [mm]: %.6f\n', obj.PleatedMembrane.Rcm);
fprintf('Rvm  [mm]: %.6f\n\n', obj.PleatedMembrane.Rvm);

fprintf('Adimensional parameters - pleated membrane\n')
fprintf('-----------------------------------------------\n')
fprintf('a1 = lv/(2*Rie):      %.2f\n', obj.PleatedMembrane.a1);
fprintf('a2 = L/lc:            %.2f\n', obj.PleatedMembrane.a2);
fprintf('a3 = Rve/Rce:         %.2f\n', obj.PleatedMembrane.a3);
fprintf('a4 = Rce/Rce_max:     %.2f\n', obj.PleatedMembrane.a4);
fprintf('a5 = Rim/Rom          %.2f\n', obj.PleatedMembrane.a5);
fprintf('a6 = lv/lc:           %.2f\n', obj.PleatedMembrane.a6);
fprintf('a7 = Rom(x=xmin)/Rib: %.2f\n\n', obj.PleatedMembrane.a7);

fprintf('Summary of design parameters - bellows \n')
fprintf('-----------------------------------------------\n')
fprintf('ns        [-]: %i\n', obj.Bellows.ns);
fprintf('L        [mm]: %.2f\n', obj.Bellows.L);
fprintf('Rib      [mm]: %.2f\n', obj.Bellows.Rib);
fprintf('Rs1      [mm]: %.2f\n', obj.Bellows.Rs1);
fprintf('Rs2      [mm]: %.2f\n', obj.Bellows.Rs2);
fprintf('lhs      [mm]: %.2f\n', obj.Bellows.lhs);
fprintf('alpha_s [deg]: %.2f\n', rad2deg(obj.Bellows.alpha_s));

if SaveSummary == 1, diary off, end

end