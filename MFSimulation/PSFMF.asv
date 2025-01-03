classdef PSFMF<handle
    properties
        offset=0;
        sigma=100; %nm
    end
    methods
        function io=intensity(obj,xpsf, xfl)
            r2=(xpsf-xfl).^2;
            rs2=sum(r2,2)/obj.sigma^2;
            io=4*exp(1)*log(2)*rs2.*exp(-4*log(2)*rs2)+obj.offset;
        end

    end
end
