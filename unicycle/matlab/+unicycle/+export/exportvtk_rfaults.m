function exportvtk_rfaults(xp,yp,zp,dim,fname,varargin)
% EXPORTVTK_RFAULTS export data to .vtp files for visualization in Paraview
%
%   exportvtk_rfaults(xp,yp,zp,dim,fname,'variable_name',data)
%
% or
%
%   exportvtk_rfaults(xp,yp,zp,fname,'variable_name1',data1,'variable_name2',data2)
%
% export a list of rectangular patches associated with data.
% 
% input:
% xp,yp,zp        - list of x,y,z coordinates of the patch corners
% dim             - number of vertices per element
% fname           - name (example.vtp) of the ASCII Paraview file
% 'variable_name' - name associated with the variable
% data            - data associated with the patch coordiantes

if 1==mod(length(varargin),2)
    error('name and associated data must come by pair')
end

fid=fopen(fname,'wt');


fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<VTKFile type=\"PolyData" version="0.1">\n');
fprintf(fid,' <PolyData>\n');

for k=1:size(xp,2)
    
    fprintf(fid,'  <Piece NumberOfPoints="%d" NumberOfPolys="1">\n',dim);
    fprintf(fid,'   <Points>\n');
    fprintf(fid,'    <DataArray type="Float32" Name="Fault Patch" NumberOfComponents="3" format="ascii">\n');
    
    % fault edge coordinates
    fprintf(fid,'%f %f %f\n',[yp(:,k) xp(:,k) -zp(:,k)]');
    
    
    fprintf(fid,'    </DataArray>\n');
    fprintf(fid,'   </Points>\n');
    fprintf(fid,'   <Polys>\n');
    fprintf(fid,'    <DataArray type="Int32" Name="connectivity" format="ascii" RangeMin="0" RangeMax="%d">\n',dim-1);
    fprintf(fid,'%d ',0:(dim-1));
    fprintf(fid,'\n');
    fprintf(fid,'    </DataArray>\n');
    fprintf(fid,'    <DataArray type="Int32" Name="offsets" format="ascii" RangeMin="%d" RangeMax="%d">\n',dim,dim);
    fprintf(fid,'%d\n',dim);
    fprintf(fid,'    </DataArray>\n');
    fprintf(fid,'   </Polys>\n');
    fprintf(fid,'   <CellData Scalars="dynamic properties">\n');
    
    for j=0:length(varargin)/2-1
        name=varargin{2*j+1};
        data=varargin{2*j+2};    
        fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n',name);
        fprintf(fid,'%f\n',data(k));
        fprintf(fid,'    </DataArray>\n');
    end
    fprintf(fid,'   </CellData>\n');
    
    fprintf(fid,'  </Piece>\n');
    
end

fprintf(fid,' </PolyData>\n');
fprintf(fid,'</VTKFile>\n');

fclose(fid);