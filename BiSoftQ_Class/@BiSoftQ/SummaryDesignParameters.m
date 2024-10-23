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
fprintf('Np     [-]: %i\n', obj.PleatedMembrane.Np);
fprintf('xtot [mm]: %.6f\n', obj.PleatedMembrane.xtot_th);
fprintf('lcp   [mm]: %.6f\n', obj.PleatedMembrane.lcp);
fprintf('lvp   [mm]: %.6f\n', obj.PleatedMembrane.lvp);
fprintf('L    [mm]: %.6f\n', obj.PleatedMembrane.L);
fprintf('Roep  [mm]: %.6f\n', obj.PleatedMembrane.Roep);
fprintf('Riep  [mm]: %.6f\n', obj.PleatedMembrane.Riep);
fprintf('Rcep  [mm]: %.6f\n', obj.PleatedMembrane.Rcep);
fprintf('Rvep  [mm]: %.6f\n', obj.PleatedMembrane.Rvep);
fprintf('Romp  [mm]: %.6f\n', obj.PleatedMembrane.Romp);
fprintf('Rimp  [mm]: %.6f\n', obj.PleatedMembrane.Rimp);
fprintf('Rcmp  [mm]: %.6f\n', obj.PleatedMembrane.Rcmp);
fprintf('Rvmp  [mm]: %.6f\n\n', obj.PleatedMembrane.Rvmp);

fprintf('Adimensional parameters - pleated membrane\n')
fprintf('-----------------------------------------------\n')
fprintf('a1 = lvp/(2*Riep):      %.2f\n', obj.PleatedMembrane.a1);
fprintf('a2 = L/lcp:            %.2f\n', obj.PleatedMembrane.a2);
fprintf('a3 = Rvep/Rcep:         %.2f\n', obj.PleatedMembrane.a3);
fprintf('a4 = Rcep/Rcep_max:     %.2f\n', obj.PleatedMembrane.a4);
fprintf('a5 = Rimp/Romp          %.2f\n', obj.PleatedMembrane.a5);
fprintf('a6 = lvp/lcp:           %.2f\n', obj.PleatedMembrane.a6);
fprintf('a7 = Romp(x=xmin)/Rib: %.2f\n\n', obj.PleatedMembrane.a7);

fprintf('Summary of design parameters - bellows \n')
fprintf('-----------------------------------------------\n')
fprintf('nb        [-]: %i\n', obj.Bellows.nb);
fprintf('L        [mm]: %.2f\n', obj.Bellows.L);
fprintf('Rib      [mm]: %.2f\n', obj.Bellows.Rib);
fprintf('Rvb      [mm]: %.2f\n', obj.Bellows.Rvb);
fprintf('Rcb      [mm]: %.2f\n', obj.Bellows.Rcb);
fprintf('lhs      [mm]: %.2f\n', obj.Bellows.lhs);
fprintf('alpha_b [deg]: %.2f\n', rad2deg(obj.Bellows.alpha_b));

if SaveSummary == 1, diary off, end

end