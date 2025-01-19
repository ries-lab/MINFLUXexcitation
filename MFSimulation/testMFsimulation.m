%MINFLUX Simulator script
L=100;
fl=fluorophore;
fl.pos=[5 3 0];
fl.brightness=100;
%initialize
psf=PSFMF_PhaseFLUX;
psf.makePSFs(L);
%initialize
loc=locMF;
loc.PSF=psf;loc.fluorophore=fl;loc.L=L;
% loc.makeorbitpattern(4,1);

%initialize
est=estMF;

numlocs=1000;
xest=zeros(numlocs,3);
phot=zeros(numlocs,1);
for k=1:numlocs
    photp=loc.patternrepeatscan;
    xest(k,:)=est.positionestimate(photp,loc.patternpos*L/2);
    phot(k)=sum(photp);
end
ff='%1.2f,';
disp(['mean(phot): ', num2str(mean(phot),ff),...
    ' std: ', num2str(std(xest),ff),...
    ' rmse: ', num2str(rmse(xest,fl.pos),ff),...
    ' bias: ', num2str(mean(xest-fl.pos),ff),...
    ' locprec: ', num2str(loc.locprecCRB(mean(phot)),ff)])
