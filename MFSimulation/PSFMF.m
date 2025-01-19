classdef PSFMF<handle
    properties

    end
    methods
        function io=intensity(obj,xpsf, xfl)
            r2=(xpsf-xfl).^2;
        end

    end
end
