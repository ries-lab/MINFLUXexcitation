classdef fluorophore<handle
    properties
        pos=[0,0,0];
        background=0;
        brightness=1; %kHz;
        blinking=false;
        ton
        toff
    end
    methods
        function Io=intensity(obj,I0)
            Io=(obj.brightness)*I0+obj.background;
        end
        function ph=photons(obj,I0)
            ph=poissrnd(obj.intensity(I0));
        end
    end
end
