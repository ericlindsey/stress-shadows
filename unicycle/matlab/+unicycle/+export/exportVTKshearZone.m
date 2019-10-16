function exportVTKshearZone(xp,yp,zp,fname,varargin)
% EXPORTVTKSHEARZONE exports data to .vtp files for visualization in Paraview
%
%   exportVTKshearZone(xp,yp,zp,dim,fname,'variable_name',data)
%
% or
%
%   exportVTKshearZone(xp,yp,zp,fname,'variable_name1',data1,'variable_name2',data2)
%
% export a list of rectangular patches associated with data.
% 
% input:
% xp,yp,zp        - list of x,y,z coordinates of the patch corners
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
    
    fprintf(fid,'  <Piece  NumberOfPoints="24" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="6">\n');
    fprintf(fid,'   <Points>\n');
    fprintf(fid,'    <DataArray type="Float32" Name="Shear Zone" NumberOfComponents="3" format="ascii">\n');
    
    % fault edge coordinates
    fprintf(fid,'%f %f %f\n',[yp(:,k) xp(:,k) -zp(:,k)]');
    
    
    fprintf(fid,'    </DataArray>\n');
    fprintf(fid,'   </Points>\n');
    fprintf(fid,'   <Polys>\n');
    fprintf(fid,'       <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="23">\n');
    fprintf(fid,'         0 1 3 2 4 6\n');
    fprintf(fid,'         7 5 8 10 11 9\n');
    fprintf(fid,'         12 13 15 14 16 18\n');
    fprintf(fid,'         19 17 20 21 23 22\n');
    fprintf(fid,'       </DataArray>\n');
    fprintf(fid,'       <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="4" RangeMax="24">\n');
    fprintf(fid,'         4 8 12 16 20 24\n');
    fprintf(fid,'       </DataArray>\n');
    fprintf(fid,'     </Polys>\n');
    
    if (~isempty(varargin))
        fprintf(fid,'   <CellData Scalars="dynamic properties">\n');
        for j=0:(length(varargin)/2-1)
            name=varargin{2*j+1};
            data=varargin{2*j+2};
            fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n',name);
            fprintf(fid,'%f\n',repmat(data(k),6,1));
            fprintf(fid,'    </DataArray>\n');
        end
        fprintf(fid,'   </CellData>\n');
    end
    
    fprintf(fid,'  </Piece>\n');
    
end

fprintf(fid,' </PolyData>\n');
fprintf(fid,'</VTKFile>\n');

fclose(fid);

    

