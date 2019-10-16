function exportxyz_rfaults(slip,x,y,z,fname,varargin)
% EXPORTXYZ_RFAULTS export data to .xyz files for processing with GMT
% tools
%
%   exportxyz_rfaults(slip,x,y,z,fname,'variable_name',data)
%
% or
%
%   exportxyz_rfaults(slip,x,y,z,fname, ...
%                     'variable_name1',data1,'variable_name2',data2)
% 
% input:
% x1,x2,x3        - list of coordinates of the patch corners (north, east, down)
% fname           - name (dir/example.flt) of the .flt Relax file
% 'variable_name' - name associated with the variable
% data            - data associated with the patch coordinates

if 1==mod(length(varargin),2)
    error('name and associated data must come by pair')
end

fid=fopen(fname,'wt');

fprintf(fid,'# export from unicycle\n');

for j=0:length(varargin)/2-1
    name=varargin{2*j+1};
    data=varargin{2*j+2};
    fprintf(fid,'# %s: %f\n',name,data);
end

for k=1:length(slip)
    fprintf(fid,'> -Z%f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n',slip(k),[x(:,k) y(:,k) z(:,k)]');
end

fclose(fid);