%o=extended(s)
%-------------
%
%Convert to extended precision.
%
%    Marcel Leutenegger © 1.1.2005
%
%  This class provides floating point computation in
%  the native floating point format of the processor.
%  The extended data type has a 64bit mantissa instead
%  of the 53bit mantissa of double data.
%
%  Due to the native execution by the processor, the
%  performance is affected only by memory throughput.
%  The current implementation minimizes memory con-
%  sumption but does not align consecutive values at
%  8 byte boundaries.
%
%  Constants and random numbers are provided by
%
%     eps(extended)     == pow2(-63)
%     pi(extended)      ~= extended(pi)
%     ones(extended)
%     rand(extended)
%     randn(extended)
%     realmax(extended)
%     realmin(extended)
%     zeros(extended)
%
%Input:
%  s     extended object, numeric or character array
%
%Output:
%  o     extended object
%