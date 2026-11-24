function position=lb_panel_position(fig,bounds,columns,index)
%LB_PANEL_POSITION Square axes in a two-row transition panel.
gapx=.026; gapy=.05; sw=(bounds(3)-(columns-1)*gapx)/columns; sh=(bounds(4)-gapy)/2;
pixels=fig.InnerPosition(3:4); side=min(sw*pixels(1),sh*pixels(2)); w=side/pixels(1); h=side/pixels(2);
column=mod(index-1,columns); row=floor((index-1)/columns);
position=[bounds(1)+column*(sw+gapx)+(sw-w)/2, ...
    bounds(2)+(1-row)*(sh+gapy)+(sh-h)/2,w,h];
end
