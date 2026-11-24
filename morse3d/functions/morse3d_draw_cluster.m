function morse3d_draw_cluster(ax,q,r0)
%MORSE3D_DRAW_CLUSTER Draw one actual 3D state as blue atoms and gray bonds.
p=reshape(q,3,[])'; p=p-mean(p,1);
hold(ax,'on');
for i=1:size(p,1)-1
    for j=i+1:size(p,1)
        if norm(p(i,:)-p(j,:))<1.25*r0
            plot3(ax,p([i j],1),p([i j],2),p([i j],3),'-','Color',[.32 .38 .45],'LineWidth',3);
        end
    end
end
[sx,sy,sz]=sphere(12);
for i=1:size(p,1)
    surf(ax,.17*sx+p(i,1),.17*sy+p(i,2),.17*sz+p(i,3), ...
        'FaceColor',[.08 .36 .88],'EdgeColor','none','FaceLighting','gouraud');
end
axis(ax,'equal'); axis(ax,'off'); view(ax,-135,28);
xlim(ax,[-1.5 1.5]);ylim(ax,[-1.5 1.5]);zlim(ax,[-1.5 1.5]);
camlight(ax,'headlight'); material(ax,'dull');
end
