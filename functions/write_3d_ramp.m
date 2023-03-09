function write_3d_ramp(patchfname, xc, yc, strike, dip, fault_width, fault_length, fault_depth, nW, nL)
    % create basic 3D ramp .seg file for Unicycle, with default parameters
    %
    % Eric Lindsey, July 2019
    
    % set fixed parameters
    v_plate = 1;
    x2 = xc -sind(strike)*fault_length/2;
    x1 = yc -cosd(strike)*fault_length/2;
    x3 = fault_depth;
    rake = 90;
    qL = 1; %parameters to increase the size of patches with distance. Uniform patches = 1
    qW = 1;
    
    patch_width = fault_width/nW; % each patch is the same size
    patch_length = fault_length/nL; % each patch is the same size
    
    % write Unicycle .seg file with one line
    fileID = fopen(patchfname,'w');
    fprintf(fileID,'%s\n', '#patch file generated automatically - 2D ramp model, constant patch size');
    fprintf(fileID,'%s\n', '# n  Vpl    x1      x2   x3   Length  Width   Strike  Dip  Rake      L0     W0    qL qW');
    fprintf(fileID,' 1  %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n', ...
        v_plate, x1, x2, x3, fault_length, fault_width, strike, dip, rake, patch_length, patch_width, qL, qW);
    fclose(fileID);
    
end
