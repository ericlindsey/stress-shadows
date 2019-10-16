function save_jointinv_model_trench0(scenario,values,fname)
    % save a 3D model, adding duplicate values along the trench x=0 
    % Eric Lindsey, May 2019

    % find indices of the patches along trenchward edge
    Ix0=find(scenario.sources{1}.geom.xc(:,1)==min(scenario.sources{1}.geom.xc(:,1)));
    % create output list: [x,y,values]
    outdata=[scenario.sources{1}.geom.xc(:,1)/1e3, scenario.sources{1}.geom.xc(:,2)/1e3, values];
    % append duplicate values along x=0
    outdata=[outdata; 0*Ix0, scenario.sources{1}.geom.xc(Ix0,2)/1e3,values(Ix0)];
    save(fname,'outdata','-ASCII'); 
    
end