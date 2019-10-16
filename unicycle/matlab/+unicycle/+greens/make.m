function make()
% function MAKE compiles C++/Fortran code EDCMP
%
%   hmmvp.kvf
%
% usage:
%   
%   hmmvp.make()
%
% note:
%   
%   for openmp, you may have to edit matlabroot/bin/mexopts.sh and change
%     
%     CC='xcrun  -sdk macosx10.9  clang'
%
%   to
% 
%     CC='/sw/bin/gcc-fsf-4.8'
%
%   and make sure 
%
%     MACOSX_DEPLOYMENT_TARGET='10.9'
%     MW_SDK_TEMP="find `xcode-select -print-path` -name MacOSX10.9.sdk"
%
%   are set to the proper values. 
%
%   Finally, edit matlabroot/bin/.matlab7rc.sh and make sure
%
%     LDPATH_PREFIX='/sw/lib/gcc4.8/lib'
%
%   to load the correct gfortran dynamic library.
%
%   mex -v is used to get the compiler in verbose mode.
%

% Configure options

p = './+unicycle/+greens';
flags = '-O CFLAGS="\$CFLAGS -fopenmp" -L../edcmp/lib -ledcmp -L/sw/lib/gcc4.8/lib -lgfortran LDFLAGS="\$LDFLAGS -fopenmp -lgomp "';

if (~isempty(dir(sprintf('%s/mexopts.bat', p))))
   mc = sprintf('mex -f %s/mexopts.bat -outdir %s %s', p, p, flags);
else
   mc = sprintf('mex -outdir %s %s', p, flags);
end

cmd=sprintf('%s ./+unicycle/+greens/calc_edcmp.cpp', mc);
fprintf('%s\n',cmd);
eval(cmd);

end
