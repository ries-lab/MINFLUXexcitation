%%

function [s] = calpitch (lambda, pitch, nth, xpix)

d = xpix * pitch; % d is the line spacing between the step structrres.
alpha = asin (nth*lambda)/d; % d(sin (alpha) = nth * lambda is the formula for blazed grating.
s = alpha;


end



