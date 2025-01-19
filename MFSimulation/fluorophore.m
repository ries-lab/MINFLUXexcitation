classdef fluorophore<handle
    properties
        pos=[0,0,0]; %nm
        background=0;
        brightness=1; %kHz;
        blinking=false;
        ton
        toff
    end
    methods
        function Io=intensity(obj,I0,brightness,background)
            if nargin<3
                brightness=obj.brightness;
            end
            if nargin<4
                background=obj.background;
            end                
            Io=(brightness)*I0+background;
        end
        function ph=photons(obj,I0,varargin)
            ph=poissrnd(obj.intensity(I0,varargin{:}));
        end
    end
end
