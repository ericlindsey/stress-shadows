function scaled_quiver(x,y,v,u,scale,plotargs)
    %deal with matlab's inability to plot multiple different quivers on the same scale.
    h = quiver(x,y,v,u,0, plotargs{:});
    hU = get(h,'UData'); 
    hV = get(h,'VData'); 
    set(h,'UData',scale*hU,'VData',scale*hV);
end
