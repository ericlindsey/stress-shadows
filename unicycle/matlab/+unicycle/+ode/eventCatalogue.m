classdef eventCatalogue < handle
    properties
        % start and end velocities
        vStart;
        vEnd;
        
        % events
        evt;
        
        % number of events
        nEvents;
        
        % is event started?
        isEvent;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = eventCatalogue(vStart,vEnd)
            % EVENT is a class representing a slip event
            %
            %   opt = ode.event(vStart,vEnd)
            %
            % where vStart and vEnd and the start and end slip velocites.
            %
            % type methods(src) and properties(src) for a list of the class
            % methods and properties.
            %
            % SEE ALSO: unicycle.
            
            obj.vStart=vStart;
            obj.vEnd=vEnd;
            
            % event catalogue
            obj.evt={};
            
            % number of events
            obj.nEvents=0;
            
            % start with no event
            obj.isEvent=false;
            
        end % constructor
        
    end % methods
    
end
