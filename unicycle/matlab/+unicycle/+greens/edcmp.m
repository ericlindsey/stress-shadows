classdef edcmp < unicycle.greens.model
    % EDCMP is a class containing static functions to compute forward
    % models between fault slip and surface observation points using the
    % semi-analytic solution of Wang et al. (2003) for a layered half space.
    methods (Static)
        function myG=G(patch,prefix,x,vecsize)
            % function G=greens.edcmp.G(src,predix,x) computes inversion kernels
            % (Green's functions) to connect fault slip to surface deformation
            % using Wang et al. (2003) analytic solution for a layered half space.
            %
            % INTERFACE:
            %
            %   G = greens.edcmp.G(src,prefix,x,vecsize)
            %
            % INPUT:
            %
            % src     geometry.patch object
            % prefix  prefix of path to the edgrn point-source Green's function
            % x       coordinates of observation point (east, north, up)
            %         *** but point should be at the surface ***
            % vecsize number of component for displacement vector
            %
            %
            % DATA LAYOUT:
            %
            %      --- strike slip ---; --- dip slip ---
            %     /e1                                   \
            %     |n1                                   |
            %     |u1                                   |
            %     |e2                                   |
            %     |n2                                   |
            %     |u2                                   |
            % G = | .                                   |
            %     | .                                   |
            %     | .                                   |
            %     | .                                   |
            %     |en                                   |
            %     |nn                                   |
            %     \un                                   /
            %
            
            assert(unicycle.greens.model.dgf==2,'greens.edcmp.G assume 2=dgf');
            
            % number of GPS stations
            D=size(x,1);
            
            % Green's function matrix
            myG=zeros(vecsize*D,unicycle.greens.model.dgf*patch.N);
            
            unicycle.greens.calc_edcmp('init',prefix);
            
            rake=[0 90];
            
            for i=1:patch.N
                % east coordinate of fault upper corner
                xg=patch.x(i,1);
                % north coordinate of fault upper corner
                yg=patch.x(i,2);
                % depth coordinate of fault upper corner
                zg=-patch.x(i,3);
                
                % length of fault path (in strike direction)
                L=patch.L(i);
                % width of fault path (in dip direction)
                W=patch.W(i);
                
                % orientation
                strike=patch.strike(i);
                dip=patch.dip(i);

                % relative distance between source and receiver
                xd=x(:,1)-xg;
                yd=x(:,2)-yg;
                
                % interlace east, north and up forward models
                for j=1:unicycle.greens.model.dgf
                    [ux,uy,uz]=unicycle.greens.calc_edcmp('d',1,0,0,zg,L,W,strike,dip,rake(j),xd,yd);
                    
                    switch(vecsize)
                        case 2
                            u=[ux,uy];
                        case 3
                            u=[ux,uy,uz];
                        otherwise
                            error('unicycle:greens:edcmp','invalid degrees of freedom (%d)\n',...
                                unicycle.greens.model.dgf);
                    end
                    for l=1:vecsize
                        % strike slip and dip slip are not interleaved
                        myG(l:vecsize:end,i+(j-1)*patch.N)=u(:,l);
                    end
                end
            end
        end % end G (Green's function)
    end % end Static methods
end % end class definition

