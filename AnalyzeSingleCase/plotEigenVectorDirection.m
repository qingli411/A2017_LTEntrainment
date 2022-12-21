function h = plotEigenVectorDirection(cmax,cmin,varargin)
% plotEigenVectorDirection(cmax,cmin,...) plots the direction of
%   the maximum or minimum eigen vector

    nArgs = length(varargin);
    if nArgs == 0
        l_wgt = 0;
    elseif nArgs == 3
        l_wgt = 1;
        lambda = varargin{1};
        wgt = varargin{2};
        l_colorbar = varargin{3};
    else
        error('plotEigenVectorDirection(cmax, cmin, [Lambda, Weight, ColorbarOn])');
    end
    
    % Vertices of the barycentric map
    xc = [-0.5, 0.5,           0];
    yc = [   0,   0, sqrt(3)*0.5];
    lc = abs(xc(2)-xc(1));
    
    % get the Cartesian coordinates
    xxmin = xc*cmin';
    yymin = yc*cmin';
    xxmax = xc*cmax';
    yymax = yc*cmax';
    
    % plot the triangle
    plot(xc, yc, 'k.');
    hold on;
    for i=1:3
        ip1 = mod(i,3)+1;
        l = line([xc(i),xc(ip1)],[yc(i),yc(ip1)]);
        l.Color = 'black';
        l.LineWidth = 2;
    end

    % grid
    nsp = 6;
    nspm1 = nsp - 1;
    r = 1.0;
    dtheta = pi/2/nsp;
    nphi = 100;
    dphi = pi/2/(nphi-1);
    cc = zeros(nphi,3);
    ctick = zeros(nspm1*3,3);
    for k = 1:3
        k1 = mod(k,3)+1;
        k2 = mod(k+1,3)+1;
        k3 = mod(k+2,3)+1;
        for i = 1:nspm1
            theta = i*dtheta;
            for j = 1:nphi
                phi = (j-1)*dphi;
                cc(j,k1) = r*cos(theta)*cos(phi);
                cc(j,k2) = r*cos(theta)*sin(phi);
                cc(j,k3) = r*sin(theta);
            end
            ctick(i+(k-1)*nspm1,:) = [cc(1,1) cc(1,2) cc(1,3)];
            cl = cc.^2;
            xl = xc*cl';
            yl = yc*cl';
            ll = plot(xl,yl,'--');
            ll.Color = 'black';
            ll.LineWidth = 1;
        end        
    end

    cticks = ctick.^2;
    xtick = xc*cticks';
    ytick = yc*cticks';

    % turn off axis
    axis off;
    % set aspect ratio to 1:1:1
    daspect([1,1,1]);

    % labels
    label = ['x';...
             'y';...
             'z'];
    clabel = cellstr(label);
    lb_pos_x = [xc(1)-0.04*lc,xc(2)+0.04*lc,xc(3)];
    lb_pos_y = [yc(1)-0.04*lc,yc(2)-0.04*lc,yc(3)+0.055*lc];
    for i=1:3
        tx = text(lb_pos_x(i),lb_pos_y(i),clabel(i));
        tx.FontSize = 16;
        tx.HorizontalAlignment = 'center';
    end
    
    % ticks and labels
    dtick = 90/nsp;
    ticks = 90-dtick:-dtick:dtick;
    ticks_str = cell(1,nspm1);
    for i = 1:nspm1
        ticks_str{i} = sprintf('$%2d^{\\circ}$',ticks(i));
    end
    % bottom
    tick_pos_x = xtick(1:nspm1);
    tick_pos_y = ytick(1:nspm1)-0.04*lc;
    for i = 1:nspm1
        tx = text(tick_pos_x(i),tick_pos_y(i),ticks_str{i});
        tx.FontSize = 12;
        tx.HorizontalAlignment = 'center';
        tx.Interpreter = 'latex';
    end
    ap1 = [tick_pos_x(end-1), tick_pos_y(end-1)-0.05*lc];
    ap2 = [tick_pos_x(2), tick_pos_y(2)-0.05*lc];
    dap = ap2-ap1;
    qv = quiver(ap1(1),ap1(2),dap(1),dap(2),0);
    qv.Color = 'black';
    qv.LineWidth = 1.5;
    tx = text(ap1(1)-0.03*lc,ap1(2),'x');
    tx.FontSize = 12;

    % right side
    tick_pos_x = xtick(nspm1+1:nspm1*2)+0.04*lc;
    tick_pos_y = ytick(nspm1+1:nspm1*2)+0.01*lc;
    for i = 1:nspm1
        tx = text(tick_pos_x(i),tick_pos_y(i),ticks_str{i});
        tx.FontSize = 12;
        tx.HorizontalAlignment = 'center';
        tx.Interpreter = 'latex';
    end
    ap1 = [tick_pos_x(end-1)+0.05*lc, tick_pos_y(end-1)+0.03*lc];
    ap2 = [tick_pos_x(2)+0.05*lc, tick_pos_y(2)+0.03*lc];
    dap = ap2-ap1;
    qv = quiver(ap1(1),ap1(2),dap(1),dap(2),0);
    qv.Color = 'black';
    qv.LineWidth = 1.5;
    tx = text(ap1(1)+0.005*lc,ap1(2)-0.02*lc,'y');
    tx.FontSize = 12;
    
    % left side
    tick_pos_x = xtick(nspm1*2+1:nspm1*3)-0.04*lc;
    tick_pos_y = ytick(nspm1*2+1:nspm1*3)+0.01*lc;
    for i = 1:nspm1
        tx = text(tick_pos_x(i),tick_pos_y(i),ticks_str{i});
        tx.FontSize = 12;
        tx.HorizontalAlignment = 'center';
        tx.Interpreter = 'latex';
    end
    ap1 = [tick_pos_x(end-1)-0.05*lc, tick_pos_y(end-1)+0.03*lc];
    ap2 = [tick_pos_x(2)-0.05*lc, tick_pos_y(2)+0.03*lc];
    dap = ap2-ap1;
    qv = quiver(ap1(1),ap1(2),dap(1),dap(2),0);
    qv.Color = 'black';
    qv.LineWidth = 1.5;
    tx = text(ap1(1)+0.005*lc,ap1(2)+0.02*lc,'z');
    tx.FontSize = 12;
     
    % plot data
    if l_wgt
        colormap(jet);
        bsize = 500;
        eta1 = lambda(:,1)-lambda(:,2);
        eta2 = 2.*(lambda(:,2)-lambda(:,3));
        inds1 = find(eta1>=eta2);
        inds2 = find(eta1<eta2);
        h = scatter(xxmax(inds1),yymax(inds1),eta1(inds1).*bsize,wgt(inds1),'o',...
                'LineWidth',1.5);
        scatter(xxmin(inds2),yymin(inds2),eta2(inds2).*bsize,wgt(inds2),'x',...
                'LineWidth',1.5);
        scatter(xc(1)+0.9*lc,0.95*yc(3),bsize,'o',...
                'LineWidth',1.5,...
                'MarkerEdgeColor','k');
        scatter(xc(1)+0.9*lc,0.85*yc(3),bsize,'x',...
                'LineWidth',1.5,...
                'MarkerEdgeColor','k');
        if l_colorbar
            cb = colorbar;
            cb.FontSize = 12;
            ylabel(cb,'z/h','FontSize',14);
        end
    else
        h = plot(xx,yy,'ko');
    end
    
end