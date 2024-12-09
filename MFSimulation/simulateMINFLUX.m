% MINFLUX simulation
%now 2D, parabolic PSF model
%times in us
%positions in nm

%TODO:
% option to add iteration if photons below a cutoff (lower than target
% photons).
%only one photon stream global. But calculate on fraction would be very
%slow, we would need to constrict search.

%test parameters
numlocs=200;
xfgt=[0.,0]; %position


%system parameters
pattern.localizationtime=1000;
pattern.patternrepeat=1;
pattern.L=50;
pattern.usecenter=false;
pattern.orbitpoints=6;
pattern.targetphot=1000;
psfpar.off=0;


% precsion in dependence on ton, toff
numpix=10;
ton=logspace(0,log10(pattern.localizationtime),numpix);
toff=logspace(0,log10(pattern.localizationtime),numpix);
bias=zeros(length(ton),length(toff),2);precision=bias;rmsenorm=bias;photons=zeros(length(ton),length(toff));
photpos=zeros(1,pattern.orbitpoints+pattern.usecenter);
for i1=1:length(ton)
    for i2=1:length(toff)
        [bias(i1,i2,:),precision(i1,i2,:),rmsenorm(i1,i2,:),photons(i1,i2),~,photall]=MFbiasprecision(xfgt,ton(i1),toff(i2),pattern,psfpar,numlocs);
        photpos=photpos+mean(photall,1);
    end
end

locprecMF=pattern.L./sqrt(8*(photons));

figure(191)
subplot(3,3,1)
imagesc(toff,ton,(rmsenorm(:,:,1)))
set(gca,'ColorScale','log')
% imagesc(toff,ton,precision(:,:,1)./locprecMF)
title(['RMSE/locprec x rel, pattern rep: ' num2str(pattern.patternrepeat)])
ylabel('ton (us)')
xlabel('toff (us)')
set(gca,{'XScale','YScale'},{'log','log'});
colorbar
axis xy


subplot(3,3,2)
imagesc(toff,ton,rmsenorm(:,:,2))
set(gca,'ColorScale','log')
% imagesc(toff,ton,precision(:,:,2)./locprecMF)
title('RMSE/locprec y rel')
ylabel('ton (us)')
xlabel('toff (us)')
set(gca,{'XScale','YScale'},{'log','log'});
colorbar
axis xy

subplot(3,3,3)
imagesc(toff,ton,bias(:,:,1))
title('bias x nm')
ylabel('ton (us)')
xlabel('toff (us)')
set(gca,{'XScale','YScale'},{'log','log'});
colorbar
axis xy

subplot(3,3,4)
imagesc(toff,ton,bias(:,:,2))
title('bias y nm')
ylabel('ton (us)')
xlabel('toff (us)')
set(gca,{'XScale','YScale'},{'log','log'});
colorbar
axis xy

subplot(3,3,5)
imagesc(toff,ton,photons)
title('photons')
ylabel('ton (us)')
xlabel('toff (us)')
set(gca,{'XScale','YScale'},{'log','log'});
colorbar
axis xy
% subplot(3,3,6)
% plot(photpos/sum(photpos))
% title('photons in each position')
% ylabel('frequency')
% xlabel('position')
photperpos=photpos/sum(photpos)
%%
numlocs=2000;
ton=100;
toff=100;
patternreps=[1:8 10:2:18 20:5:40 50:10:100];
bias=zeros(length(patternreps),2);precision=bias;rmsenorm=bias; photons=zeros(length(patternreps),1);
xest=zeros(length(patternreps),numlocs,2);
photpos=zeros(1,pattern.orbitpoints+pattern.usecenter);
for k=1:length(patternreps)
    pattern.patternrepeat=patternreps(k);
    [bias(k,:),precision(k,:), rmsenorm(k,:),photons(k),xest(k,:,:),photall]=MFbiasprecision(xfgt,ton,toff,pattern,psfpar,numlocs);
    photpos=photpos+mean(photall,1);
end
locprecMF=pattern.L./sqrt(8*(photons(:)));
figure(191)
subplot(3,3,7)
% plot(patternreps,precision./locprecMF)
plot(patternreps,rmsenorm)
title(['RMSE/locprec, ton=' num2str(ton) ', toff=' num2str(toff)])
ylabel('RMSE/locprec rel')
xlabel('pattern repeats')
ax=gca;
ax.YLim(1)=1;

subplot(3,3,8)
plot(patternreps,bias)
title('bias')
ylabel('bias nm')
xlabel('pattern repeats')


subplot(3,3,9)
figure(77)
plot(xest([1 end],:,1)')
xlabel('data point')
ylabel('x position nm')

photperpos=photpos/sum(photpos)
%%
pattern2=pattern;
pattern2.patternrepeat=100;
[bias,precision, rmsenorm,photons,xest,photall]=MFbiasprecision(xfgt,100,100,pattern2,psfpar,1000);
ff=3;
disp(['x: ' num2str(bias(1),ff) ' ± ' num2str(precision(1),ff) ' nm, y: ' num2str(bias(2),ff) ' ± ' num2str(precision(2),ff) ' nm, locprec ' num2str(mean(locprecMF(:)),ff) 'nm'])


%%

function [bias,precision,rmsenorm,photons,xest,photall]=MFbiasprecision(xfgt,ton,toff,pattern,psfpar,numlocs)
for k=numlocs:-1:1
    [xest(k,:),phot(k),photall(k,:)]=simulateMFphotons(xfgt,ton,toff,pattern,psfpar);
end
dx=xest-xfgt;
precision=std(dx,[],1,'omitmissing');

locprecMF=pattern.L./sqrt(8*(phot(:)));
rmsenorm=sqrt(mean(dx.^2./locprecMF.^2,'omitmissing'));

% rmse=sqrt(mean(dx.^2,'omitmissing'));
bias=mean(dx,1,'omitmissing');
photons=mean(phot,'omitmissing');
% sum(isnan(bias))
% mean(photall,1)
end

function [xest,numphot,photp]=simulateMFphotons(xfgt,ton,toff,pattern,psfpar)
%calculate derived parameters
totalpoints=pattern.orbitpoints+pattern.usecenter;
time_pattern_point=pattern.localizationtime/pattern.patternrepeat/totalpoints;
pattern.xpos=makepattern(pattern);
io=PSF([0 pattern.L/2],psfpar);
itot=io*pattern.orbitpoints*pattern.patternrepeat;
fraction=ton/(ton+toff);
% fraction=1;
brightness=pattern.targetphot/itot/fraction;

maxsimultime=totalpoints*time_pattern_point*pattern.patternrepeat;

tonoff=blinkingtrace(ton,toff,maxsimultime);
time_pattern_start=rand*tonoff(end,end)/2; %choose random point to avoid starting always with the on state
%but is this right? we coarsly localize when in on state. Should be option.
%repeat if too few photons?
inten=zeros(1,totalpoints);
photpint=zeros(1,totalpoints);

for pr=1:pattern.patternrepeat
    iPSF=PSF(-(xfgt-pattern.xpos(:,:)),psfpar);
    for p=1:totalpoints
        fraction_on=calculate_onfraction(tonoff,time_pattern_start,time_pattern_point);
        % inten(p)=PSF(-(xfgt-pattern.xpos(p,:)),psfpar)*fraction_on;
        photpint(p)=photpint(p)+iPSF(p)*fraction_on*brightness;
        time_pattern_start=time_pattern_start+time_pattern_point;
    end
end
photp=poissrnd(photpint);
xest=estimate_position(photp,pattern,psfpar);
numphot=sum(photp);
end



function tonoff=blinkingtrace(ton,toff,maxsimultime)
numpoints=max(10,ceil(maxsimultime/(ton+toff)*4));
% numpoints=ceil(maxsimultime/(ton+toff)*2);
tont=exprnd(ton,numpoints,1);
tofft=exprnd(toff,numpoints,1);
tsum=tont+tofft;
timestart=cumsum([0;tsum(1:end-1)]);
tonoff=horzcat(tont,tofft,timestart);
if timestart(end)<maxsimultime*2 %too short simulation: add another one
    tonoffx=blinkingtrace(ton,toff,maxsimultime);
    tonoff=vertcat(tonoff,tonoffx);
    tsum=sum(tonoff(:,1:2),2);
    tonoff(:,3)=cumsum([0;tsum(1:end-1)]);
end
end

function xpattern=makepattern(pattern)
dphi=2*pi/pattern.orbitpoints;
phi=(0:dphi:2*pi-dphi)';
x=cos(phi)*pattern.L/2;
y=sin(phi)*pattern.L/2;
xpattern=horzcat(x,y);
if pattern.usecenter
    xpattern(end+1,:)=[0,0];
end
end

function fraction=calculate_onfraction(tonoff,t0,dwelltime)
tbstart=tonoff(:,3);
ton=tonoff(:,1);
toff=tonoff(:,2);
i0=find(t0<tbstart,1,"first")-1;
i1=find(t0+dwelltime<tbstart,1,"first")-1;
if isempty(i1)
    disp('blinking trace too short')
    tonoff, t0, dwelltime
end
ontime=dwelltime-sum(toff(i0:i1-1));
t0rel=t0-tbstart(i0);
if t0rel>ton(i0) %correct for first off
    ontime=ontime+(t0rel-ton(i0));
end
t1rel=t0+dwelltime-tbstart(i1);
if t1rel>ton(i1) %correct for first off
    ontime=ontime-(t1rel-ton(i1));
end
fraction=max(ontime/dwelltime,0);
end


function testfraction
z=ones(5,1);
tonoff=horzcat(z,z,(0:2:8)');
fraction=calculate_onfraction(tonoff,0,3)
end


function intensity=PSF(xrel,psfpar)
intensity=(sum(xrel.^2,2))+psfpar.off;
end


function xest=estimate_position(photp,pattern,psfpar)
pi=photp/sum(photp);
% eq 2.63
xest=-sum(pi'.*pattern.xpos);
end