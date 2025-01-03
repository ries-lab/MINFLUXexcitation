classdef estMF<handle
    properties
        PSF

    end
    methods
        function xest=positionestimate(obj,photonsi,patternpos)
            pi=photonsi/sum(photonsi);
            % eq 2.63
            xest=-sum(pi'.*patternpos);
        end

    end
end
