function plot_field_slices(x,y,z,vel,bfl,xpos)
% Plots slice (cross planes) in a 3D plot

    arg = permute(vel,[3 1 2]);
    %arg = flip(arg,1);

    [X,Z,Y] = meshgrid(x,z,y);

    numslices = length(xpos);

    for i = 1:numslices

        slice(X,Z,Y,arg,x(xpos(i)),[],[],'spline'); hold on;
        set(findobj(gca,'Type','Surface'),'EdgeColor','none');

        colormap(redblue(512)); 
        if i ==1
            cax = caxis;
            caxis([-1 1]*max(abs(cax)));
        else
            caxis([-1 1]*max(abs(cax)));
        end

        ypos = find(bfl(xpos(i),:)>=0.9999,1,'first');

        plot3([x(xpos(i)) x(xpos(i))],[-25 25],[y(ypos) y(ypos)],...
            '--k','LineWidth',1.0);

    end

    daspect([1 1 0.5]);
    yticks([z(1) 0 ceil(z(end))]);
    axis([x(1) x(800) z(1) ceil(z(end)) 0 10]);
    set(gca,'Ydir','reverse');

    box on;
    xlabel('x');ylabel('z');zlabel('y');
end