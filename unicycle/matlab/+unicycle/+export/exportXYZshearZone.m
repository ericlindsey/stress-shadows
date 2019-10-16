function exportXYZshearZone(value,x,y,z,fname)
% EXPORTXYZSHEARZONE export data to .xyz files for processing with GMT
% tools
%
%   exportXYZshearZone(x,y,z,fname,'variable_name',data)
%
% or
%
%   exportXYZshearZone(x,y,z,fname, ...
%                     'variable_name1',data1,'variable_name2',data2)
% 
% input:
% value           - for the -Z option in GMT .xyz file
% x1,x2,x3        - list of coordinates of the patch corners (north, east, down)
% fname           - name (dir/example.flt) of the .flt Relax file

fid=fopen(fname,'wt');

fprintf(fid,'# export from unicycle\n');

for k=1:length(value)
    for j=1:6
        fprintf(fid,'> -Z%f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n', ...
            value(k),[x((1+4*(j-1)):4*j,k),y((1+4*(j-1)):4*j,k),z((1+4*(j-1)):4*j,k)]');
    end
end

fclose(fid);