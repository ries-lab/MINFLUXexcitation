%% This function adds a phase shift to the entire area of the input image.
% Code written by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.
% This is meant to be used with MATLAB library, "Electromagnetic field in the focus region of a microscope objective", 
% based on the original publication, Leutenegger et al., "Fast focus field calculations", published in Optics Express 14 (2006) 

function E = piston(sys,E,r,t,p)

E = E.*exp((sys.pistonphi)*i);

end

