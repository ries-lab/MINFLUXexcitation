%MINFLUX Simulator script

fl=fluorophore;
fl.brightness=100;
%initialize
psf=PSFMF;
%initialize
loc=locMF;
loc.PSF=psf;
loc.fluorophore=fl;
loc.L=100;
loc.makeorbitpattern(4,1);
%initialize
est=estMF;

numlocs=1000;
xest=zeros(numlocs,3);
phot=zeros(numlocs,1);
for k=1:numlocs
    photp=loc.patternrepeatscan;
    xest(k,:)=est.positionestimate(photp,loc.patternpos);
    phot(k)=sum(photp);
end
std(xest)
loc.locprecCRB(mean(phot))
mean(phot)