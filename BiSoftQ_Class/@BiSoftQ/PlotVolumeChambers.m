function PlotVolumeChambers(obj)

% create title string
mytitle_line1 = strcat("$N_p$ = ", num2str(obj.PleatedMembrane.Np)," , ", ... 
                       "$L/l_{c,p}$ = ", num2str(obj.PleatedMembrane.a2), ", ", ...
                       "$l_{v,p}/l_{c,p}$ = ", num2str(obj.PleatedMembrane.a6),", ", ...
                       "$R_{im,p}/R_{om,p}$ = ", num2str(round(obj.PleatedMembrane.a5,2)), ", ");

mytitle_line2 = strcat("$l_{v,p}/(2\,R_{ie,p})$ = ", num2str(round(obj.PleatedMembrane.a1,2)),", ", ...
                       "$R_{i,b}$ = ", num2str(round(obj.PleatedMembrane.Rib,1))," mm, ", ...
                       "$R_{ie,p}$ = ", num2str(round(obj.PleatedMembrane.Riep,1)), " mm ");

% plot volume
fig(1) = figure('Name','Volume');
hold on
plot(obj.x_vec,obj.VolumeChamber1_vec*1e-3)
plot(obj.x_vec,obj.VolumeChamber2_vec*1e-3)
hold off
xlabel('x (mm)')
ylabel('Volume (cm$^3$)')
grid on
box on
legend('Chamber 1','Chamber 2','interpreter','latex','Location','best')
title({mytitle_line1; mytitle_line2})

end