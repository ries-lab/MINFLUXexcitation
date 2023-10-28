%Normally distributed random numbers.
%
%	Marcel Leutenegger © 2.5.2005
%
%
%This function uses the MATLAB generators for normally
%and uniformly distributed random numbers. It provides
%random numbers with 64 bit mantissa. The most signifi-
%cant 53 bits are normally distributed. The lasting 11
%bits are uniformly distributed.
%
%Hence, the normal distribution is interpolated piece-
%wise linear.
%
%The function can repeat a series of random numbers if
%*both* MATLAB random number generators are initialised
%to a user-defined value *and* the series of the output
%array sizes are identical.
%