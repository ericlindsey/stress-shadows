classdef evt < handle
    properties
        % start and end times
        tStart;
        tEnd;
        
        % start and end indices
        iStart;
        iEnd;
        
        % solution vector at start and end times
        yStart;
        yEnd;
        
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = evt(tStart,tEnd,iStart,iEnd,yStart,yEnd)
            % EVENT is a class representing a slip event
            %
            %   opt = ode.event(tStart,tEnd,iStart,iEnd,yStart,yEnd)
            %
            % where xStart and xEnd and the start and end x attributes.
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle.
            
            obj.tStart=tStart;
            obj.tEnd=tEnd;
            
            obj.iStart=iStart;
            obj.iEnd=iEnd;
            
            obj.yStart=yStart;
            obj.yEnd=yEnd;
            
        end % constructor
        
    end % methods
    
end
