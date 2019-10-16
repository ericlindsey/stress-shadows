function plot_faults(faults,dim,xlim,ylim,varargin)
% function PLOT_FAULTS(faults,dim,xlim,ylim [,color])
%

if 4<nargin
    color=varargin{1};
else
    color=[0.5 0.5 0.5];
end

hold on

[sfl,~]=size(dim);
for i=1:sfl
    if i==1
        S1=1;
    else
        S1=S2+1;
    end
    S2=S1+dim(i)-1;
    ab=faults(S1:S2,1);
    or=faults(S1:S2,2);
    
    pos=ab>xlim(1) & ab<xlim(2) & or>ylim(1) & or<ylim(2);
    if 0~=numel(pos)
        line(ab(pos),or(pos),'Color',color,'LineStyle','-','LineWidth',0.5); 
    end
end
