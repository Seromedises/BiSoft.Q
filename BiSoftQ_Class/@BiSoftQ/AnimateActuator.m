function AnimateActuator(obj,options)

arguments
    obj BiSoftQ 
    options.NPoints (1,1) double = 100
    options.SaveVideo (1,1) logical = false
    options.VideoName (1,:) string = "BiSoftQ_Animation"
    options.SeparateVideos (1,1) logical = false
end

NPoints = options.NPoints;
SaveVideo = options.SaveVideo;
VideoName = options.VideoName;
SeparateVideos = options.SeparateVideos;

L = obj.PleatedMembrane.L;
N = obj.PleatedMembrane.N;
Roe = obj.PleatedMembrane.Roe;
Rib = obj.PleatedMembrane.Rib;

colors = [0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980];

% initial configuration
% pleated membrane
[~,~,~,~,~,X_membr,Y_membr,Z_membr] = obj.PleatedMembrane.PleatedMembraneSurface(...
    obj.PleatedMembrane.rc0,obj.PleatedMembrane.rc,obj.PleatedMembrane.rs1, ...
    obj.PleatedMembrane.x0s1,obj.PleatedMembrane.rs2,obj.PleatedMembrane.x0s2,obj.PleatedMembrane.vmax);

% bellows
[X_bell,Y_bell,Z_bell] = obj.Bellows.BellowsSurface(obj.Bellows.alpha_s,obj.Bellows.h, ...
    obj.Bellows.rcs1,obj.Bellows.Rs1,obj.Bellows.Rt1,obj.Bellows.h1, ...
    obj.Bellows.rcs2,obj.Bellows.Rs2,obj.Bellows.Rt2,0);

% preallocate matrices
X_membr_frame = NaN([size(X_membr),obj.NPoints,2]);
Y_membr_frame = NaN([size(Y_membr),obj.NPoints,2]);
Z_membr_frame = NaN([size(Z_membr),obj.NPoints,2]);

X_bell_frame = NaN([size(X_bell),obj.NPoints,2]);
Y_bell_frame = NaN([size(Y_bell),obj.NPoints,2]);
Z_bell_frame = NaN([size(Z_bell),obj.NPoints,2]);

for jj=1:2
    % initial configuration frame
    X_membr_frame(:,:,1,jj) = X_membr;
    Y_membr_frame(:,:,1,jj) = Y_membr;
    Z_membr_frame(:,:,1,jj) = Z_membr;

    X_bell_frame(:,:,1,jj) = X_bell;
    Y_bell_frame(:,:,1,jj) = Y_bell;
    Z_bell_frame(:,:,1,jj) = Z_bell;
end

if isempty(obj.F_pull_vec)
    % compute force characteristic
    obj.ForceCharacteristic("NPoints",NPoints,"method","method_2","approx_method","",...
        "filetype","mfile","ElasticForceBool",true,"NumWorkers",0);
end

index = (length(obj.F_pull_vec)+1)/2;
F_pull_max = 3*obj.F_pull_vec(index);

x_vec = [obj.x_vec(1:index); obj.x_vec(index); obj.x_vec(index+1:end)];
x_vec = reshape(x_vec,index,2);
x_vec(:,1) = flip(x_vec(:,1));

F_pull_vec = [obj.F_pull_vec(1:index); obj.F_pull_vec(index); obj.F_pull_vec(index+1:end)];
F_pull_vec = reshape(F_pull_vec,index,2);
F_pull_vec(:,1) = flip(F_pull_vec(:,1));

F_push_vec = [obj.F_push_vec(1:index); obj.F_push_vec(index); obj.F_push_vec(index+1:end)];
F_push_vec = reshape(F_push_vec,index,2);
F_push_vec(:,1) = flip(F_push_vec(:,1));

idxs = zeros(2,1);
for ii=1:2
    if isempty(find(isnan(x_vec(:,ii)),1,'first'))
        idxs(ii) = length(x_vec(:,ii))+1;
    else
        idxs(ii) = find(isnan(x_vec(:,ii)),1,'first');
    end
end

C = struct2cell(obj.PleatedMembrane.DeformedConfigsParameters);
[rc0_vec,rc_vec,~,~,rs1_vec,x0s1_vec,rs2_vec,x0s2_vec,vmax_vec] = C{:};

C = struct2cell(obj.Bellows.DeformedConfigsParameters);
[alpha_s_vec,h_vec,rcs1_vec,Rs1_vec,Rt1_vec,h1_vec,rcs2_vec,Rs2_vec,Rt2_vec] = C{:};

% deformed configuration
iter_max = size(x_vec,1);
parfor jj=1:2
    
    for ii=2:iter_max
        
        if ii==idxs(jj)
            break 
        else

            [~,~,~,~,~,X_membr,Y_membr,Z_membr] = obj.PleatedMembrane.PleatedMembraneSurface(rc0_vec(ii,jj), ...
                rc_vec(ii,jj),rs1_vec(ii,jj),x0s1_vec(ii,jj),rs2_vec(ii,jj),x0s2_vec(ii,jj),vmax_vec(ii,jj));
          
            [X_bell,Y_bell,Z_bell] = obj.Bellows.BellowsSurface(alpha_s_vec(ii,jj),h_vec(ii,jj), ...
                rcs1_vec(ii,jj),Rs1_vec(ii,jj),Rt1_vec(ii,jj),h1_vec(ii,jj),rcs2_vec(ii,jj), ...
                Rs2_vec(ii,jj),Rt2_vec(ii,jj),x_vec(ii,jj));
            
            X_membr_frame(:,:,ii,jj) = X_membr;
            Y_membr_frame(:,:,ii,jj) = Y_membr;
            Z_membr_frame(:,:,ii,jj) = Z_membr-x_vec(ii,jj)/2;
            
            X_bell_frame(:,:,ii,jj) = X_bell;
            Y_bell_frame(:,:,ii,jj) = Y_bell;
            Z_bell_frame(:,:,ii,jj) = Z_bell;

        end

    end

end

if strcmp(obj.type,"A")
    XMax = max(abs(max(X_membr_frame,[],'all')),abs(min(X_membr_frame,[],'all')));
    YMax = max(abs(max(Y_membr_frame,[],'all')),abs(min(Y_membr_frame,[],'all')));
    ZMax = max(abs(max(Z_membr_frame,[],'all')),abs(min(Z_membr_frame,[],'all')));
    FaceColorBell = [0.5 0.5 0.5];
    FaceAlphaBell = 1;
    FaceColorMembr = [0.90 0.90 0.90];
    FaceAlphaMembr = 0.9;
elseif strcmp(obj.type,"B")
    XMax = max(abs(max(X_bell_frame,[],'all')),abs(min(X_bell_frame,[],'all')));
    YMax = max(abs(max(Y_bell_frame,[],'all')),abs(min(Y_bell_frame,[],'all')));
    ZMax = max(abs(max(Z_bell_frame,[],'all')),abs(min(Z_bell_frame,[],'all')));
    FaceColorBell = [0.90 0.90 0.90];
    FaceAlphaBell = 0.9;
    FaceColorMembr = [0.5 0.5 0.5];
    FaceAlphaMembr = 1;
end

mytitle_line = ['BiSoft.Q type ',obj.type,' - isometric view'];
 
if SaveVideo==1
    fig(1) = figure('Name','Phase 1','units','pixels','position',[0 0 1920 1080],'ToolBar','none','Color','white');
else
    fig(1) = figure('Name','Phase 1','units','pixels','position',[0 0 1920 1080],'ToolBar','none');
end

ii = 1;
% from x=0 to x=xmin, pulling
for pp=1:obj.NPoints+1

    if isnan(X_membr_frame(1,1,pp,1))
        idx_end_1 = pp-1;
        break
    end

    X_membr_def = squeeze(X_membr_frame(:,:,pp,1));
    Y_membr_def = squeeze(Y_membr_frame(:,:,pp,1));
    Z_membr_def = squeeze(Z_membr_frame(:,:,pp,1));

    X_bell_def = squeeze(X_bell_frame(:,:,pp,1));
    Y_bell_def = squeeze(Y_bell_frame(:,:,pp,1));
    Z_bell_def = squeeze(Z_bell_frame(:,:,pp,1));

    % 3d plot of the parametric surface
    subplot(2,2,[1,3],'Parent',fig(1))
    hold off
    % bellows
    surf(X_bell_def,Y_bell_def,Z_bell_def,'EdgeColor','none','FaceColor',FaceColorBell,'FaceAlpha',FaceAlphaBell)
    hold on
    plot3(X_bell_def(1,:),Y_bell_def(1,:),Z_bell_def(1,:),'k-','LineWidth',1.5)
    for jj=2:size(Z_bell_def,1)/50   
        plot3(X_bell_def((jj-1)*50,:),Y_bell_def((jj-1)*50,:),Z_bell_def((jj-1)*50,:),'k-','LineWidth',1.5)
    end
    plot3(X_bell_def(end,:),Y_bell_def(end,:),Z_bell_def(end,:),'k-','LineWidth',1.5)
    % pleated membrane
    surf(X_membr_def,Y_membr_def,Z_membr_def,'EdgeColor','none','FaceColor',FaceColorMembr,'FaceAlpha',FaceAlphaMembr);
    for iii=1:2*N
        idx = (iii-1)*25+1;
        plot3(X_membr_def(idx,:),Y_membr_def(idx,:),Z_membr_def(idx,:),'k','LineWidth',1.5)
    end
    plot3(X_membr_def(:,1),Y_membr_def(:,1),Z_membr_def(:,1),'k','LineWidth',1.5)
    % plot3(X_membr_def(:,25),Y_membr_def(:,25),Z_membr_def(:,25),'k','LineWidth',1.5)
    plot3(X_membr_def(:,end),Y_membr_def(:,end),Z_membr_def(:,end),'k','LineWidth',1.5)
    hold off
    axis equal
    grid off
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    xlim([-XMax XMax])
    ylim([-YMax YMax])
    zlim([-ZMax ZMax])
    % view(0,90)
    view(35,30)
    title(mytitle_line)
    
    subplot(2,2,2,'Parent',fig(1))
    hold off
    rho_c = sqrt(X_membr_frame(1,:,pp,1).^2+Y_membr_frame(1,:,pp,1).^2);
    rho_v = sqrt(X_membr_frame(25,:,pp,1).^2+Y_membr_frame(25,:,pp,1).^2);

    if strcmp(obj.type,"A")

        plot(1.3*[-max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all') max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')],...
        [-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.3*max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')+1,-L/2+2,'$x=0$','interpreter','latex')
        
        hold on

        fill([X_bell_def(1,1) rho_c X_bell_def(end,1) flip(X_bell_def(:,1))'],[Z_bell_def(1,1) Z_membr_frame(1,:,pp,1) Z_bell_def(end,1) flip(Z_bell_def(:,1))'],'c','FaceAlpha',0.5)
        fill([-X_bell_def(1,1) -rho_c -X_bell_def(end,1) -flip(X_bell_def(:,1))'],[Z_bell_def(1,1) Z_membr_frame(1,:,pp,1) Z_bell_def(end,1) flip(Z_bell_def(:,1))'],'c','FaceAlpha',0.5)
   
    elseif strcmp(obj.type,"B")

        plot(1.20*[-XMax XMax],[-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.20*XMax+1,-L/2+1,'$x=0$','interpreter','latex')
        
        hold on

        fill([-rho_c(1) -rho_c -rho_c(end) flip(rho_c) -rho_c(1)],...
        [Z_membr_frame(1,1,pp,1) Z_membr_frame(1,:,pp,1) Z_membr_frame(1,end,pp,1) flip(Z_membr_frame(1,:,pp,1)) Z_membr_frame(1,1,pp,1)],...
        'c','FaceAlpha',0.5)

    end
    
    hold on

    H1 = plot(rho_c,Z_membr_frame(1,:,pp,1),'b');
    plot(-rho_c,Z_membr_frame(1,:,pp,1),'b')

    H2 = plot(rho_v,Z_membr_frame(25,:,pp,1),'r');            
    plot(-rho_v,Z_membr_frame(25,:,pp,1),'r')

    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,end,pp,1) Z_membr_frame(1,end,pp,1)],'b')
    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,1,pp,1) Z_membr_frame(1,1,pp,1)],'b')

    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,end,pp,1) Z_membr_frame(1,end,pp,1)],'r')
    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,1,pp,1) Z_membr_frame(1,1,pp,1)],'r')

    R_sec = sqrt(X_bell_def(:,1).^2+Y_bell_def(:,1).^2);
    [R_sec_mirror,Z_sec_mirror] = obj.PleatedMembrane.MirrorPoints(pi/2,R_sec,Z_bell_def(:,1));
    R_sec = [R_sec;R_sec_mirror;R_sec(1)];
    Z_sec = [Z_bell_def(:,1);Z_sec_mirror;Z_bell_def(1,1)];
    plot(R_sec,Z_sec,'k-')

    if strcmp(obj.type,"A")
        n = 24;
        DD = Roe/n;
        plot([-Roe -Roe+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Roe+2*(iii-1)*DD -Roe+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    elseif strcmp(obj.type,"B")
        n = 24;
        DD = Rib/n;
        plot([-Rib -Rib+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Rib+2*(iii-1)*DD -Rib+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    end

    axis equal
    hold off
    legend([H1 H2],'crest','valley','Location','northeastoutside','interpreter','latex')
    xlim(1.3*[-XMax XMax])
    ylim([-ZMax ZMax])
    box on
    title('Longitudinal sections - Chamber 1 active')
    
    subplot(2,2,4,'Parent',fig(1))
    hold off
    H1 = plot(x_vec(1:pp,1),F_pull_vec(1:pp,1),'-','Color',colors(1,:));
    hold on
    H2 = plot(NaN,NaN,'-','Color',colors(2,:));
    xlim([min(x_vec,[],'all') max(x_vec,[],'all')])
    ylim([min(2*F_push_vec,[],'all') F_pull_max])
    hold off
    xlabel('x (mm)')
    ylabel('F (N)')
    if strcmp(obj.type,"A")
        legend([H1 H2],'Chamber 1 active','Chamber 2 active','Location','northwest','interpreter','latex')
    elseif strcmp(obj.type,"B")
        legend([H1 H2],'Chamber 1 active','Chamber 1+2 active','Location','northwest','interpreter','latex')
    end
    box on
    title(['Force characteristic (p = ',num2str(obj.p_gauge_MPa*10),' bar)'])
  
    drawnow
    M(ii,1) = getframe(gcf);
    ii = ii+1;

    idx_end_1 = pp;
   
end

if SaveVideo==1
    fig(2) = figure('Name','Phase 2','units','pixels','position',[0 0 1920 1080],'ToolBar','none','Color','white');
else
    fig(2) = figure('Name','Phase 2','units','pixels','position',[0 0 1920 1080],'ToolBar','none');
end

ii = 1;
% from x=xmin to x=0, pushing
for pp=idx_end_1:-1:1

    if isnan(X_membr_frame(1,1,pp,1))
        continue
    end

    X_membr_def = squeeze(X_membr_frame(:,:,pp,1));
    Y_membr_def = squeeze(Y_membr_frame(:,:,pp,1));
    Z_membr_def = squeeze(Z_membr_frame(:,:,pp,1));

    X_bell_def = squeeze(X_bell_frame(:,:,pp,1));
    Y_bell_def = squeeze(Y_bell_frame(:,:,pp,1));
    Z_bell_def = squeeze(Z_bell_frame(:,:,pp,1));

    % 3d plot of the parametric surface
    subplot(2,2,[1,3],'Parent',fig(2))
    % bellows
    surf(X_bell_def,Y_bell_def,Z_bell_def,'EdgeColor','none','FaceColor',FaceColorBell,'FaceAlpha',FaceAlphaBell)
    hold on
    plot3(X_bell_def(1,:),Y_bell_def(1,:),Z_bell_def(1,:),'k-','LineWidth',1.5)
    for jj=2:size(Z_bell_def,1)/50   
        plot3(X_bell_def((jj-1)*50,:),Y_bell_def((jj-1)*50,:),Z_bell_def((jj-1)*50,:),'k-','LineWidth',1.5)
    end
    plot3(X_bell_def(end,:),Y_bell_def(end,:),Z_bell_def(end,:),'k-','LineWidth',1.5)
    % pleated membrane
    surf(X_membr_def,Y_membr_def,Z_membr_def,'EdgeColor','none','FaceColor',FaceColorMembr,'FaceAlpha',FaceAlphaMembr);
    for iii=1:2*N
        idx = (iii-1)*25+1;
        plot3(X_membr_def(idx,:),Y_membr_def(idx,:),Z_membr_def(idx,:),'k','LineWidth',1.5)
    end
    plot3(X_membr_def(:,1),Y_membr_def(:,1),Z_membr_def(:,1),'k','LineWidth',1.5)
    % plot3(X_membr_def(:,25),Y_membr_def(:,25),Z_membr_def(:,25),'k','LineWidth',1.5)
    plot3(X_membr_def(:,end),Y_membr_def(:,end),Z_membr_def(:,end),'k','LineWidth',1.5)
    hold off
    axis equal
    grid off
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    xlim([-XMax XMax])
    ylim([-YMax YMax])
    zlim([-ZMax ZMax])
    % view(0,90)
    view(35,30)
    title(mytitle_line)
    
    subplot(2,2,2,'Parent',fig(2))
    hold off
    rho_c = sqrt(X_membr_frame(1,:,pp,1).^2+Y_membr_frame(1,:,pp,1).^2);
    rho_v = sqrt(X_membr_frame(25,:,pp,1).^2+Y_membr_frame(25,:,pp,1).^2);
    
    R_sec = sqrt(X_bell_def(:,1).^2+Y_bell_def(:,1).^2);
    [R_sec_mirror,Z_sec_mirror] = obj.PleatedMembrane.MirrorPoints(pi/2,R_sec,Z_bell_def(:,1));
    R_sec = [R_sec;R_sec_mirror;R_sec(1)];
    Z_sec = [Z_bell_def(:,1);Z_sec_mirror;Z_bell_def(1,1)];
    
    if strcmp(obj.type,"A")

        plot(1.3*[-max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all') max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')],...
        [-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.3*max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')+1,-L/2+2,'$x=0$','interpreter','latex')
        
        hold on

        fill([X_bell_def(1,1) X_bell_def(:,1)' X_bell_def(end,1) -X_bell_def(end,1) -flip(X_bell_def(:,1))' -X_bell_def(1,1)],...
        [Z_bell_def(1,1) Z_bell_def(:,1)' Z_bell_def(end,1) Z_bell_def(end,1) flip(Z_bell_def(:,1))' Z_bell_def(1,1)],'c','FaceAlpha',0.5)
    
    elseif strcmp(obj.type,"B")

        plot(1.20*[-XMax XMax],[-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.20*XMax+1,-L/2+1,'$x=0$','interpreter','latex')
        
        hold on

        fill(R_sec,Z_sec,'c','FaceAlpha',0.5)

    end
    
    hold on

    plot(R_sec,Z_sec,'k-')

    H1 = plot(rho_c,Z_membr_frame(1,:,pp,1),'b');
    plot(-rho_c,Z_membr_frame(1,:,pp,1),'b')
   
    H2 = plot(rho_v,Z_membr_frame(25,:,pp,1),'r');            
    plot(-rho_v,Z_membr_frame(25,:,pp,1),'r')

    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,end,pp,1) Z_membr_frame(1,end,pp,1)],'b')
    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,1,pp,1) Z_membr_frame(1,1,pp,1)],'b')

    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,end,pp,1) Z_membr_frame(1,end,pp,1)],'r')
    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,1,pp,1) Z_membr_frame(1,1,pp,1)],'r')

    if strcmp(obj.type,"A")
        n = 24;
        DD = Roe/n;
        plot([-Roe -Roe+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Roe+2*(iii-1)*DD -Roe+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    elseif strcmp(obj.type,"B")
        n = 24;
        DD = Rib/n;
        plot([-Rib -Rib+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Rib+2*(iii-1)*DD -Rib+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    end

    axis equal
    hold off
    legend([H1 H2],'crest','valley','Location','northeastoutside','interpreter','latex')
    xlim(1.3*[-XMax XMax])
    ylim([-ZMax ZMax])
    box on
    if strcmp(obj.type,"A")
        title('Longitudinal sections - Chamber 2 active')
    elseif strcmp(obj.type,"B")
        title('Longitudinal sections - Chamber 1+2 active')
    end

    subplot(2,2,4,'Parent',fig(2))
    hold on
    H2 = plot(x_vec(pp:idx_end_1,1),F_push_vec(pp:idx_end_1,1),'-','Color',colors(2,:));
    H1 = plot(x_vec(:,1),F_pull_vec(:,1),'-','Color',colors(1,:));
    xlim([min(x_vec,[],'all') max(x_vec,[],'all')])
    ylim([min(2*F_push_vec,[],'all') F_pull_max])
    hold off
    xlabel('x (mm)')
    ylabel('F (N)')
    if strcmp(obj.type,"A")
        legend([H1 H2],'Chamber 1 active','Chamber 2 active','Location','northwest','interpreter','latex')
    elseif strcmp(obj.type,"B")
        legend([H1 H2],'Chamber 1 active','Chamber 1+2 active','Location','northwest','interpreter','latex')
    end
    box on
    title(['Force characteristic (p = ',num2str(obj.p_gauge_MPa*10),' bar)'])
  
    drawnow
    M(ii,2) = getframe(gcf);
    ii = ii+1;
    
end

if SaveVideo==1
    fig(3) = figure('Name','Phase 3','units','pixels','position',[0 0 1920 1080],'ToolBar','none','Color','white');
else
    fig(3) = figure('Name','Phase 3','units','pixels','position',[0 0 1920 1080],'ToolBar','none');
end

ii = 1;
% from x=0 to x=xmax, pushing
for pp=1:obj.NPoints+1

    if isnan(X_membr_frame(1,1,pp,2))
        break
    end

    X_membr_def = squeeze(X_membr_frame(:,:,pp,2));
    Y_membr_def = squeeze(Y_membr_frame(:,:,pp,2));
    Z_membr_def = squeeze(Z_membr_frame(:,:,pp,2));

    X_bell_def = squeeze(X_bell_frame(:,:,pp,2));
    Y_bell_def = squeeze(Y_bell_frame(:,:,pp,2));
    Z_bell_def = squeeze(Z_bell_frame(:,:,pp,2));

    % 3d plot of the parametric surface
    subplot(2,2,[1,3],'Parent',fig(3))
    hold off
    % bellows
    surf(X_bell_def,Y_bell_def,Z_bell_def,'EdgeColor','none','FaceColor',FaceColorBell,'FaceAlpha',FaceAlphaBell)
    hold on
    plot3(X_bell_def(1,:),Y_bell_def(1,:),Z_bell_def(1,:),'k-','LineWidth',1.5)
    for jj=2:size(Z_bell_def,1)/50   
        plot3(X_bell_def((jj-1)*50,:),Y_bell_def((jj-1)*50,:),Z_bell_def((jj-1)*50,:),'k-','LineWidth',1.5)
    end
    plot3(X_bell_def(end,:),Y_bell_def(end,:),Z_bell_def(end,:),'k-','LineWidth',1.5)
    % pleated membrane
    surf(X_membr_def,Y_membr_def,Z_membr_def,'EdgeColor','none','FaceColor',FaceColorMembr,'FaceAlpha',FaceAlphaMembr);
    for iii=1:2*N
        idx = (iii-1)*25+1;
        plot3(X_membr_def(idx,:),Y_membr_def(idx,:),Z_membr_def(idx,:),'k','LineWidth',1.5)
    end
    plot3(X_membr_def(:,1),Y_membr_def(:,1),Z_membr_def(:,1),'k','LineWidth',1.5)
    % plot3(X_membr_def(:,25),Y_membr_def(:,25),Z_membr_def(:,25),'k','LineWidth',1.5)
    plot3(X_membr_def(:,end),Y_membr_def(:,end),Z_membr_def(:,end),'k','LineWidth',1.5) 
    axis equal
    grid off
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    xlim([-XMax XMax])
    ylim([-YMax YMax])
    zlim([-ZMax ZMax])
    % view(0,90)
    view(35,30)
    title(mytitle_line)
    
    subplot(2,2,2,'Parent',fig(3))
    hold off
    rho_c = sqrt(X_membr_frame(1,:,pp,2).^2+Y_membr_frame(1,:,pp,2).^2);
    rho_v = sqrt(X_membr_frame(25,:,pp,2).^2+Y_membr_frame(25,:,pp,2).^2);

    R_sec = sqrt(X_bell_def(:,1).^2+Y_bell_def(:,1).^2);
    [R_sec_mirror,Z_sec_mirror] = obj.PleatedMembrane.MirrorPoints(pi/2,R_sec,Z_bell_def(:,1));
    R_sec = [R_sec;R_sec_mirror;R_sec(1)];
    Z_sec = [Z_bell_def(:,1);Z_sec_mirror;Z_bell_def(1,1)];

    if strcmp(obj.type,"A")
        
        plot(1.3*[-max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all') max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')],...
        [-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.3*max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')+1,-L/2+2,'$x=0$','interpreter','latex')
        
        hold on

        fill([X_bell_def(1,1) X_bell_def(:,1)' X_bell_def(end,1) -X_bell_def(end,1) -flip(X_bell_def(:,1))' -X_bell_def(1,1)],...
        [Z_bell_def(1,1) Z_bell_def(:,1)' Z_bell_def(end,1) Z_bell_def(end,1) flip(Z_bell_def(:,1))' Z_bell_def(1,1)],'c','FaceAlpha',0.5)
    
    elseif strcmp(obj.type,"B")

        plot(1.20*[-XMax XMax],[-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.20*XMax+1,-L/2+1,'$x=0$','interpreter','latex')
        
        hold on

        fill(R_sec,Z_sec,'c','FaceAlpha',0.5)

    end

    hold on

    plot(R_sec,Z_sec,'k-')

    H1 = plot(rho_c,Z_membr_frame(1,:,pp,2),'b');
    plot(-rho_c,Z_membr_frame(1,:,pp,2),'b')

    H2 = plot(rho_v,Z_membr_frame(25,:,pp,2),'r');            
    plot(-rho_v,Z_membr_frame(25,:,pp,2),'r')

    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,end,pp,2) Z_membr_frame(1,end,pp,2)],'b')
    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,1,pp,2) Z_membr_frame(1,1,pp,2)],'b')

    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,end,pp,2) Z_membr_frame(1,end,pp,2)],'r')
    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,1,pp,2) Z_membr_frame(1,1,pp,2)],'r')

    if strcmp(obj.type,"A")
        n = 24;
        DD = Roe/n;
        plot([-Roe -Roe+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Roe+2*(iii-1)*DD -Roe+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    elseif strcmp(obj.type,"B")
        n = 24;
        DD = Rib/n;
        plot([-Rib -Rib+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Rib+2*(iii-1)*DD -Rib+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    end

    axis equal
    hold off
    legend([H1 H2],'crest','valley','Location','northeastoutside','interpreter','latex')
    xlim(1.3*[-XMax XMax])
    ylim([-ZMax ZMax])
    box on
    if strcmp(obj.type,"A")
        title('Longitudinal sections - Chamber 2 active')
    elseif strcmp(obj.type,"B")
        title('Longitudinal sections - Chamber 1+2 active')
    end

    subplot(2,2,4,'Parent',fig(3))
    hold on
    H2 = plot(x_vec(1:pp,2),F_push_vec(1:pp,2),'-','Color',colors(2,:));
    plot(x_vec(:,1),F_push_vec(:,1),'-','Color',colors(2,:))
    H1 = plot(x_vec(:,1),F_pull_vec(:,1),'-','Color',colors(1,:));
    xlim([min(x_vec,[],'all') max(x_vec,[],'all')])
    ylim([min(2*F_push_vec,[],'all') F_pull_max])
    hold off
    xlabel('x (mm)')
    ylabel('F (N)')
    if strcmp(obj.type,"A")
        legend([H1 H2],'Chamber 1 active','Chamber 2 active','Location','northwest','interpreter','latex')
    elseif strcmp(obj.type,"B")
        legend([H1 H2],'Chamber 1 active','Chamber 1+2 active','Location','northwest','interpreter','latex')
    end
    box on
    title(['Force characteristic (p = ',num2str(obj.p_gauge_MPa*10),' bar)'])
  
    drawnow
    M(ii,3) = getframe(gcf);
    ii = ii+1;

end

if SaveVideo==1
    fig(4) = figure('Name','Phase 4','units','pixels','position',[0 0 1920 1080],'ToolBar','none','Color','white');
else
    fig(4) = figure('Name','Phase 4','units','pixels','position',[0 0 1920 1080],'ToolBar','none');
end

ii = 1;
% from x=xmax to x=0, pulling
for pp=obj.NPoints+1:-1:1

    if isnan(X_membr_frame(1,1,pp,2))
        continue
    end

    X_membr_def = squeeze(X_membr_frame(:,:,pp,2));
    Y_membr_def = squeeze(Y_membr_frame(:,:,pp,2));
    Z_membr_def = squeeze(Z_membr_frame(:,:,pp,2));

    X_bell_def = squeeze(X_bell_frame(:,:,pp,2));
    Y_bell_def = squeeze(Y_bell_frame(:,:,pp,2));
    Z_bell_def = squeeze(Z_bell_frame(:,:,pp,2));

    % 3d plot of the parametric surface
    subplot(2,2,[1,3],'Parent',fig(4))
    hold off
    % bellows
    surf(X_bell_def,Y_bell_def,Z_bell_def,'EdgeColor','none','FaceColor',FaceColorBell,'FaceAlpha',FaceAlphaBell)
    hold on
    plot3(X_bell_def(1,:),Y_bell_def(1,:),Z_bell_def(1,:),'k-','LineWidth',1.5)
    for jj=2:size(Z_bell_def,1)/50   
        plot3(X_bell_def((jj-1)*50,:),Y_bell_def((jj-1)*50,:),Z_bell_def((jj-1)*50,:),'k-','LineWidth',1.5)
    end
    plot3(X_bell_def(end,:),Y_bell_def(end,:),Z_bell_def(end,:),'k-','LineWidth',1.5)
    % pleated membrane
    surf(X_membr_def,Y_membr_def,Z_membr_def,'EdgeColor','none','FaceColor',FaceColorMembr,'FaceAlpha',FaceAlphaMembr);
    for iii=1:2*N
        idx = (iii-1)*25+1;
        plot3(X_membr_def(idx,:),Y_membr_def(idx,:),Z_membr_def(idx,:),'k','LineWidth',1.5)
    end
    plot3(X_membr_def(:,1),Y_membr_def(:,1),Z_membr_def(:,1),'k','LineWidth',1.5)
    % plot3(X_membr_def(:,25),Y_membr_def(:,25),Z_membr_def(:,25),'k','LineWidth',1.5)
    plot3(X_membr_def(:,end),Y_membr_def(:,end),Z_membr_def(:,end),'k','LineWidth',1.5) 
    axis equal
    grid off
    xlabel('X (mm)')
    ylabel('Y (mm)')
    zlabel('Z (mm)')
    xlim([-XMax XMax])
    ylim([-YMax YMax])
    zlim([-ZMax ZMax])
    % view(0,90)
    view(35,30)
    title(mytitle_line)
    
    subplot(2,2,2,'Parent',fig(4))
    hold off
    rho_c = sqrt(X_membr_frame(1,:,pp,2).^2+Y_membr_frame(1,:,pp,2).^2);
    rho_v = sqrt(X_membr_frame(25,:,pp,2).^2+Y_membr_frame(25,:,pp,2).^2);
    
    if strcmp(obj.type,"A")

        plot(1.3*[-max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all') max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')],...
        [-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.3*max(sqrt(X_membr_frame(1,:,:,:).^2+Y_membr_frame(1,:,:,:).^2),[],'all')+1,-L/2+2,'$x=0$','interpreter','latex')
        
        hold on
    
        fill([X_bell_def(1,1) rho_c X_bell_def(end,1) flip(X_bell_def(:,1))'],[Z_bell_def(1,1) Z_membr_frame(1,:,pp,2) Z_bell_def(end,1) flip(Z_bell_def(:,1))'],'c','FaceAlpha',0.5)
        fill([-X_bell_def(1,1) -rho_c -X_bell_def(end,1) -flip(X_bell_def(:,1))'],[Z_bell_def(1,1) Z_membr_frame(1,:,pp,2) Z_bell_def(end,1) flip(Z_bell_def(:,1))'],'c','FaceAlpha',0.5)
        
    elseif strcmp(obj.type,"B")

        plot(1.20*[-XMax XMax],[-L/2 -L/2],'k--','Linewidth',0.75)
        text(-1.20*XMax+1,-L/2+1,'$x=0$','interpreter','latex')
        
        hold on

        fill([-rho_c(1) -rho_c -rho_c(end) flip(rho_c) -rho_c(1)],...
        [Z_membr_frame(1,1,pp,2) Z_membr_frame(1,:,pp,2) Z_membr_frame(1,end,pp,2) flip(Z_membr_frame(1,:,pp,2)) Z_membr_frame(1,1,pp,2)],...
        'c','FaceAlpha',0.5)

    end
    
    hold on

    H1 = plot(rho_c,Z_membr_frame(1,:,pp,2),'b');
    plot(-rho_c,Z_membr_frame(1,:,pp,2),'b')

    H2 = plot(rho_v,Z_membr_frame(25,:,pp,2),'r');            
    plot(-rho_v,Z_membr_frame(25,:,pp,2),'r')

    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,end,pp,2) Z_membr_frame(1,end,pp,2)],'b')
    plot([-rho_c(1) rho_c(1)],[Z_membr_frame(1,1,pp,2) Z_membr_frame(1,1,pp,2)],'b')

    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,end,pp,2) Z_membr_frame(1,end,pp,2)],'r')
    plot([-rho_v(1) rho_v(1)],[Z_membr_frame(1,1,pp,2) Z_membr_frame(1,1,pp,2)],'r')

    R_sec = sqrt(X_bell_def(:,1).^2+Y_bell_def(:,1).^2);
    [R_sec_mirror,Z_sec_mirror] = obj.PleatedMembrane.MirrorPoints(pi/2,R_sec,Z_bell_def(:,1));
    R_sec = [R_sec;R_sec_mirror;R_sec(1)];
    Z_sec = [Z_bell_def(:,1);Z_sec_mirror;Z_bell_def(1,1)];
    plot(R_sec,Z_sec,'k-')

    if strcmp(obj.type,"A")
        n = 24;
        DD = Roe/n;
        plot([-Roe -Roe+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Roe+2*(iii-1)*DD -Roe+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    elseif strcmp(obj.type,"B")
        n = 24;
        DD = Rib/n;
        plot([-Rib -Rib+2*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        hold on
        for iii=2:n+1
            plot([-Rib+2*(iii-1)*DD -Rib+2*iii*DD],[L/2 L/2+2*DD],'k-','Linewidth',0.75)
        end
    end

    axis equal
    hold off
    legend([H1 H2],'crest','valley','Location','northeastoutside','interpreter','latex')
    xlim(1.3*[-XMax XMax])
    ylim([-ZMax ZMax])
    box on
    title('Longitudinal sections - Chamber 1 active')

    subplot(2,2,4,'Parent',fig(4))
    hold on
    plot(x_vec(pp:obj.NPoints+1,2),F_pull_vec(pp:obj.NPoints+1,2),'-','Color',colors(1,:));
    plot(x_vec,F_push_vec,'-','Color',colors(2,:));
    H1 = plot(NaN,NaN,'-','Color',colors(1,:));
    H2 = plot(NaN,NaN,'-','Color',colors(2,:));
    plot(x_vec(:,1),F_pull_vec(:,1),'-','Color',colors(1,:))
    xlim([min(x_vec,[],'all') max(x_vec,[],'all')])
    ylim([min(2*F_push_vec,[],'all') F_pull_max])
    hold off
    xlabel('x (mm)')
    ylabel('F (N)')
    if strcmp(obj.type,"A")
        legend([H1 H2],'Chamber 1 active','Chamber 2 active','Location','northwest','interpreter','latex')
    elseif strcmp(obj.type,"B")
         legend([H1 H2],'Chamber 1 active','Chamber 1+2 active','Location','northwest','interpreter','latex')
    end
    box on
    title(['Force characteristic (p = ',num2str(obj.p_gauge_MPa*10),' bar)'])
  
    drawnow
    M(ii,4) = getframe(gcf);
    ii = ii+1;        
end

if SaveVideo==1 && SeparateVideos==0
    FileName = VideoName;
    v = VideoWriter(FileName,'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 5;
    open(v)
    for jj=1:4
        for ii=1:size(M,1)
            if isempty(M(ii,jj).cdata)
                continue
            end
            writeVideo(v,M(ii,jj));
        end
    end
    close(v)
end

if SaveVideo==1 && SeparateVideos==1
    for jj=1:4
        FileName = strcat(VideoName,'_',num2str(jj));
        v = VideoWriter(FileName,'MPEG-4');
        v.Quality = 100;
        fr = abs(obj.xmin_eff/obj.xmax_eff);
        if fr>=1
            if jj==3 || jj==4 
                v.FrameRate = 30*fr;
            else
                v.FrameRate = 30;
            end
        else
            if jj==1 || jj==2 
                v.FrameRate = 30/fr;
            else
                v.FrameRate = 30;
            end 
        end
        open(v)
        for ii=1:size(M,1)
            if isempty(M(ii,jj).cdata)
                continue
            end
            writeVideo(v,M(ii,jj));
        end
        close(v)
    end
end


end