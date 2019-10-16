classdef trianglePassiveReceiver < unicycle.geometry.triangle
    properties
        % static friction
        mu0;
        
        % shear stress in strike direction due to strike slip
        Kss;
        % shear stress in strike direction due to dip slip
        Kds;
        % shear stress in dip direction due to strike slip
        Ksd;
        % shear stress in dip direction due to dip slip
        Kdd;
        
        % stress components: S.xx, S.xy, S.xz, S.yy, S.yz, S.zz, S.tn
        % stress components due to strike slip
        SS;
        % stress components due to dip slip
        DS;
        
        % modeled shear stress in strike direction due to rcv and source
        sr,ss;
        % modeled shear stress in dip direction due to receiver and source
        dr,ds;
        
        % modeled tension due to receiver and source
        tr,ts;
        % modeled pressure due to receiver and source
        pr,ps;
        
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = trianglePassiveReceiver()
            % TRIANGLEPASSIVERECEIVER is a class representing the geometry
            % and physicalproperties of receiver fault patches, including 
            % position, orientation, dimension and friction properties.
            %
            %   src = geometry.trianglePassiveReceiver()
            %
            % creates a instance of fault triangles for a passive receiver 
            % fault for stress change calculations.
            %
            % SEE ALSO: unicycle
            
            if (0==nargin)
                return
            end
        end
        
        function [] = simulation(obj,evl)
            % method SIMULATION builds forward models of stress if
            % simulation exists
            
            obj.sr=zeros(obj.N,numel(evl.t));
            obj.dr=zeros(obj.N,numel(evl.t));
            obj.tr=zeros(obj.N,numel(evl.t));
            obj.pr=zeros(obj.N,numel(evl.t));
            
            for k=1:length(evl.t)
                % receiver shear stress in strike direction model
                obj.sr(:,k)=obj.Kss*evl.y(1:evl.dgf:end,k)+...
                            obj.Kds*evl.y(2:evl.dgf:end,k);
                % receiver pressure model due to strike slip and dip slip
                obj.dr(:,k)=obj.Ksd*evl.y(1:evl.dgf:end,k)+...
                            obj.Kdd*evl.y(2:evl.dgf:end,k);
                        
                % receiver tension model due to strike slip and dip slip
                obj.tr(:,k)=obj.SS.Ktn*evl.y(1:evl.dgf:end,k)+...
                            obj.DS.Ktn*evl.y(2:evl.dgf:end,k);
                % receiver pressure model due to strike slip and dip slip
                obj.pr(:,k)=-(obj.SS.Kxx+obj.SS.Kyy+obj.SS.Kzz)*evl.y(1:evl.dgf:end,k)+...
                            -(obj.DS.Kxx+obj.DS.Kyy+obj.DS.Kzz)*evl.y(2:evl.dgf:end,k);
            end
        end % simulation
        
    end % methods
    
end