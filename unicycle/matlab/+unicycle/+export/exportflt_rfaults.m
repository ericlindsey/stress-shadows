function exportflt_rfaults(slip,x1,x2,x3,len,width,strike,dip,rake,fname,varargin)
% EXPORTFLT_RFAULTS export data to .flt files for processing with Relax
% tools
%
%   exportflt_rfaults(slip,x1,x2,x3,len,width,strike,dip,rake,fname, ...
%                     'variable_name',data)
%
% or
%
%   exportflt_rfaults(slip,x1,x2,x3,len,width,strike,dip,rake,fname, ...
%                     'variable_name1',data1,'variable_name2',data2)
% 
% input:
% slip            - list of slip values
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
fprintf(fid,'# n x1 x2 x3 length width strike dip rake\n');
fprintf(fid,'%d %f %f %f %f %f %f %f %f %f\n',[(1:length(slip))',slip,x1,x2,x3,len,width,strike,dip,rake]');

fclose(fid);