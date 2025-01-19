classdef locMF<handle
    properties
        posc=[0 0 0];
        L=100; %nanometers
        patternpos=[0 0 0];
        fluorophore
        PSF
        dwelltime=1;%milliseconds
        patternrepeat=1;
    end
    methods
        function phot=patternrepeatscan(obj)
            phot=obj.orbitscan;
            for k=2:obj.patternrepeat
                phot=phot+obj.orbitscan;
            end
        end
        function phot=orbitscan(obj)
            totalpoints=size(obj.patternpos,1);
            photint=zeros(1,totalpoints);
            for k=1:totalpoints
                photint(k)=photint(k)+obj.fluorophore.intensity(obj.PSF.intensity([1,1],(obj.patternpos(k,:)*obj.L/2-obj.posc),obj.fluorophore.pos))*obj.dwelltime;
            end
            phot=poissrnd(photint);
        end
        function makeorbitpattern(obj,orbitpoints,usecenter)
            dphi=2*pi/orbitpoints;
            phi=(0:dphi:2*pi-dphi)';
            x=cos(phi);
            y=sin(phi);
            xpattern=horzcat(x,y,0*phi);
            if usecenter
                xpattern(end+1,:)=[0,0,0];
            end
            obj.patternpos=xpattern;
        end
        function lp=locprecCRB(obj,photons)
            lp=obj.L./sqrt(8*(photons));
        end
    end
end
