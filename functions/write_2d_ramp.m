function write_2d_ramp(patchfname, dip, fault_width, npatch)
    % create basic 2D ramp .seg file for Unicycle, with default parameters
    %
    % Eric Lindsey, June 2019
    
    % set fixed parameters
    v_plate = 1;
    x1 = -5e8;
    x2 = 0;
    x3 = 0;
    fault_length = 1e9; % we are using 3D greens functions (Okada) to simulate a 2d fault. So just make it long, 1e6 km
    strike = 0;
    rake = 90;
    patch_length = 1e9; % same as fault length, so there is only one patch in this direction
    qL = 1; %parameters to increase the size of patches with distance. Uniform patches = 1
    qW = 1;
    patch_width = fault_width/npatch; % each patch is the same size
    
    % write Unicycle .seg file with one line
    fileID = fopen(patchfname,'w');
    fprintf(fileID,'%s\n', '#patch file generated automatically - 2D ramp model, constant patch size');
    fprintf(fileID,'%s\n', '# n  Vpl    x1      x2   x3   Length  Width   Strike  Dip  Rake      L0     W0    qL qW');
    fprintf(fileID,' 1  %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n', ...
        v_plate, x1, x2, x3, fault_length, fault_width, strike, dip, rake, patch_length, patch_width, qL, qW);
    fclose(fileID);
    
end
