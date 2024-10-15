function PlotForceCharacteristic(obj)

% create title string
mytitle_line1 = strcat("N = ", num2str(obj.PleatedMembrane.N)," , ", ... 
                       "$L/l_c$ = ", num2str(obj.PleatedMembrane.a2), ", ", ...
                       "$l_v/l_c$ = ", num2str(obj.PleatedMembrane.a6),", ", ...
                       "$R_{i,m}/R_{o,m}$ = ", num2str(round(obj.PleatedMembrane.a5,2)), ", ");

mytitle_line2 = strcat("$l_v/(2\,R_{i,e})$ = ", num2str(round(obj.PleatedMembrane.a1,2)),", ", ...
                       "$R_{i,b}$ = ", num2str(round(obj.PleatedMembrane.Rib,1))," mm, ", ...
                       "$R_{i,e}$ = ", num2str(round(obj.PleatedMembrane.Rie,1)), " mm ");

% plot force characteristic for a set of geometric parameters
% F_pull_max = 3*obj.F_pull_vec((length(obj.F_pull_vec)+1)/2);
F_pull_max = max(obj.F_pull_vec);
fig(1) = figure('Name','Force characteristic - complete model');
hold on
P1 = plot(obj.x_vec,obj.F_pull_vec,'--');
color_plot = P1.Color;
plot(obj.x_vec,obj.F_push_vec,'-','Color',color_plot)
xlabel('x (mm)')
ylabel('F (N)')
ylim([-inf F_pull_max])
if strcmp(obj.type,"A")
    legend('Chamber 1 active','Chamber 2 active','interpreter','latex','Location','northwest')
elseif strcmp(obj.type,"B")
    legend('Chamber 1 active','Chamber 1+2 active','interpreter','latex','Location','northwest')
end
grid on
box on
hold off
title({mytitle_line1; mytitle_line2})

end