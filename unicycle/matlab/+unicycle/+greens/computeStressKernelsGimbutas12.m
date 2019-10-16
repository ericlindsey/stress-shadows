function [Kss,Ksd,Kds,Kdd]=computeStressKernelsGimbutas12(src,rcv,G,nu)
% COMPUTESTRESSKERNELSGIMBUTAS12 evalutes the stress kernels for triangle
% dislocations using the method of Gimbutas et al. (2012) using the HMMVP
% interface of Bradley.
%
% INPUT:
% src   - source triangle faults (IGNORED)
% rcv   - receiver triangle faults
% G     - rigidity
% nu    - Poisson's ratio
%
% ** Currently on the receiver's self interactions are implemented. **
%
% Sylvain Barbot, 03/30/2015
% Earth Observatory of Singapore
%
% SEE ALSO: unicycle, hmmvp

import hmmvp.*

h=unihmmvp('buildGreensFunctionsGimbutas12',src,rcv,G,nu,'./hmmvp/TGF_',1e-6,12);

% receiver stress interactions with HMMVP v1.3
[ss.id,ss.nnz]=hmmvp('init',h.Kss.hm_filename);
[sd.id,ss.nnz]=hmmvp('init',h.Ksd.hm_filename);
[ds.id,ss.nnz]=hmmvp('init',h.Kds.hm_filename);
[dd.id,ss.nnz]=hmmvp('init',h.Kdd.hm_filename);


% extract full matrices
is=1:rcv.N;
js=1:rcv.N;

Kss=hmmvp('extract',ss.id,is,js);
Ksd=hmmvp('extract',sd.id,is,js);
Kds=hmmvp('extract',ds.id,is,js);
Kdd=hmmvp('extract',dd.id,is,js);

% free the link to the hierarchy matrices
hmmvp('cleanup',ss.id);
hmmvp('cleanup',sd.id);
hmmvp('cleanup',ds.id);
hmmvp('cleanup',dd.id);

end
