function [KK,Ksn,Kdn] = unicycle_stress_kernel(source,slipComponents,stressKernelFolder)
    
    % Unicycle stress kernel computation
    % returns either the strike/dip component matrix only, or the full 2D
    % stress (traction) kernel matrix.
    % slipComponents should be 1 (strike), 2 (dip) or [1,2] (both).
    
    % note that it now also returns the Ksn and Kdn matrix separately,
    % for normal stress change calculations.
    
    % check for the existence of a saved HDF5 file containing the matrix.
    % the filename must be unique to this source object -
    % simplest way is just to use a hash of the data in the source object. 
    fname = [stressKernelFolder '/stresskernel_' DataHash(source) '.h5'];
    
    
    if (exist(fname, 'file'))
        disp(['Loading existing stress kernel matrices from file: ' fname])
        Kss=h5read(fname,'/Kss');
        Kds=h5read(fname,'/Kds');
        Ksd=h5read(fname,'/Ksd');
        Kdd=h5read(fname,'/Kdd');
        Ksn=h5read(fname,'/Ksn');
        Kdn=h5read(fname,'/Kdn');
    else
        % create the matrix if no file is found
    
        % compute the 2D stress kernel first (all 4 shear traction components)
        disp('Precomputed stress kernel matrices not found for this fault geometry. Calculating a new set...')
        tic
        [Kss,Kds,Ksd,Kdd,Ksn,Kdn]=source.tractionKernels(source);
        toc
        
        % save data to the HDF5 file so it doesn't have to be computed again
        disp(['Saving stress kernel matrices to file: ' fname])
        h5create(fname,'/Kss',size(Kss));
        h5create(fname,'/Kds',size(Kss));
        h5create(fname,'/Ksd',size(Kss));
        h5create(fname,'/Kdd',size(Kss));
        h5create(fname,'/Ksn',size(Kss));
        h5create(fname,'/Kdn',size(Kss));
        h5write(fname,'/Kss',Kss);
        h5write(fname,'/Kds',Kds);
        h5write(fname,'/Ksd',Ksd);
        h5write(fname,'/Kdd',Kdd);
        h5write(fname,'/Ksn',Ksn);
        h5write(fname,'/Kdn',Kdn);
    end
    
    % create the final matrix
    if isequal(slipComponents, [1,2])
        % Jointinv expects a single K matrix.
        % We simply place the components next to each other, in keeping
        % with the structure of the modelVector and designMatrix.
        KK = [Kss Kds; Ksd Kdd];
        
    elseif slipComponents==1
        KK = Kss;

    elseif slipComponents==2
        KK = Kdd;

    else
        error('Error: slipComponents value not recognized! Must be 1 (strike), 2 (dip) or [1,2] (both).')
    end
    
end