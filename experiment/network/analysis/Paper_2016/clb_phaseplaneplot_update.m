function clb_phaseplaneplot_update(obj,ed,vh,ph,mV,x,y,bl1)
rv = round(obj.Value);
[obj.Value, rv]
set(vh,'XData',rv*50,'YData',mV(rv*50));
set(ph,'XData',x(rv*50),'YData',y(rv*50));
set(bl1, 'String', num2str(rv*50));
drawnow;
end
