%% This function introduces a phase shift at the center of the input image.
% Code written by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.
% This is meant to be used with MATLAB library, "Electromagnetic field in the focus region of a microscope objective", 
% based on the original publication, Leutenegger et al., "Fast focus field calculations", published in Optics Express 14 (2006) 

function E = pishift(sys,E,r,t,p)

radshift = 1/sqrt(2)+0.03;
idx = r<=radshift;       % The area for pi shift.
idx = 2*idx-1;
phase = idx.*exp((pi)*i);   
if isfield(sys, 'delshift')
    phase(r<=radshift) = phase(r<=radshift).*exp(sys.delshift*i);
end
phase = repmat(phase,size(E,1),1);
E = E.*phase;
end

