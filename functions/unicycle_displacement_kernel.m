function G_3d = unicycle_displacement_kernel(geom,coords,slipComponents,kernelFolder) 
    % Unicycle Green's functions calculation
    % output G matrix with 3 components of displacement interleaved by point
    % With model parameters strike and dip (not interleaved), or only one
    % Points must be on the surface (z=0); coordinates are assumed to be in meters.
    %
    % First searches for an existing displacement kernel under the folder:
    % [kernelFolder '/displacementkernel_' DataHash({geom,coords}) '.h5']
    %
    % Runs the following unicycle computation only when no existing kernel file is found:
    % [Gs,Gd]=geom.earthModel.displacementKernels(geom, [coords(:,1),coords(:,2),0*coords(:,1)], 3);
    %
    % Eric Lindsey, 2018

    
    % check for the existence of a saved HDF5 file containing the matrix.
    % the filename must be unique to this combination of source object and GPS data -
    % simplest way is just to use a hash of the data.
    fname = [kernelFolder '/displacementkernel_' DataHash({geom,coords}) '.h5'];
    
    if (exist(fname, 'file'))
        disp(['Loading existing displacement kernel matrices from file: ' fname])
        Gs=h5read(fname,'/Gs');
        Gd=h5read(fname,'/Gd');
    else
        % create the matrix if no file is found
        disp('Precomputed displacement kernel matrices not found for this combination of')
        disp('fault geometry and observation coordinates. Calculating a new set...')
        % compute the 2D displacement kernel first (two slip components)
        tic
        [Gs,Gd]=geom.earthModel.displacementKernels(geom, [coords(:,1),coords(:,2),0*coords(:,1)], 3);
        toc
        
        % save data to the HDF5 file so it doesn't have to be computed again
        disp(['Saving displacement kernel matrices to file: ' fname])
        h5create(fname,'/Gs',size(Gs));
        h5create(fname,'/Gd',size(Gd));
        h5write(fname,'/Gs',Gs);
        h5write(fname,'/Gd',Gd);
    end
    
    % now keep only the requested components
    if isequal(slipComponents, [1,2])
        % Jointinv expects a single G matrix. We put the Strike and Dip matrices side-by-side.
        G_3d = [Gs Gd];
    elseif slipComponents==1
        % strike only
        G_3d = Gs;
    elseif slipComponents==2
        % dip only
        G_3d = Gd;
    else
        error('Error: slipComponents value not recognized! Must be 1 (strike), 2 (dip) or [1,2] (both).')
    end
    
end

% old version: abandoned because it was too confusing
% Interleave Gs, Gd by column, to keep the  (strike, dip)
% components of each patch next to each other.
% so the output G is now:
% [patch1strike_pt1E patch1dip_pt1E ... patchNdip_pt1E ...]
% [patch1strike_pt1N patch1dip_pt1N ... patchNdip_pt1N ...]
% [patch1strike_pt1U patch1dip_pt1U ... patchNdip_pt1U ...]
%G_3d = reshape([Gs; Gd], size(Gs,1), []);
