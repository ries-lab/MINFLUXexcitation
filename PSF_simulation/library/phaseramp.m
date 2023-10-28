%Phase ramp.
%
% sys    System data
%
%
% E      Electric field [V/m]
% r      Radial position [Da/2]
% t      Incidence angle [rad]
% p      Polar angle [rad]
%
function E = phaseramp(sys,E,r,t,p)

% % %%
% % %Added with Giuseppe in order to introduce a shift of the doughnut hole at
% % %the back aperture of the objective lens.
% % 
% [xx,yy] = pol2cart(p,r);
% %xx = xx+0;
% %yy = yy - 2/sys.Da*sys.OY; % Introducing a spatial offset to the phase in y direction. Use this when Ei is { 'gauss', 'offset','phaseramp',  'zernike', 'circular'};
% %yy = yy - 2/sys.Da*sys.oy(1); % Use this when Ei is { 'gauss', 'phaseramp', 'offset', 'zernike', 'circular'};
% [p,r] = cart2pol(xx,yy); 
% % %till above line.


%%
% creation of 2pi helical ramp of vortex
%e = exp(1i.*(p+pi/4));
e = exp(1i.*(p));
E=E.*repmat(e,size(E,1),1);     % Extending a matrix dimension by repeating 2D matrix 3 times in the new first dimension.



%% %Plotting the Intensity and Phase
% figure
% fig2=(squeeze(E(1,:,:)));
% imagesc(abs(fig2))
% colorbar
% axis square
% figure
% imagesc(angle(fig2))
% colorbar
% axis square


E;

%}

%%


end
 




