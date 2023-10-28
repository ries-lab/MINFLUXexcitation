
				     MATLAB toolbox
		________________________________________________________

				Extended precision class
		________________________________________________________


** Contents

	1. Introduction
	2. Requirements
	3. Installation
	4. Comments
	5. Copyright
	6. Warranty
	7. History
	8. Download
	9. Trademarks

** Publisher

	Marcel Leutenegger		marcel.leutenegger@a3.epfl.ch
	MPI for Biophysical Chemistry
	NanoBiophotonics
	Am Fassberg 11			Phone:	+49 551 201 2516
	DE-37077 Göttingen



1. Introduction

	This class provides floating point calculation in the native floating point
	format of the processor. The extended class stores data with 64bit mantissa
	instead of 53bit in double.

	The class is intended for calculations depending critically on the mantissa
	length. Summing a large with a small value generally leads to a truncation
	error of the small value to fit into the mantissa. The extended precision
	relaxes this constraint, useful to verify the accuracy
	of a numerical result.


2. Requirements

	• An Intel Pentium II compatible computer or newer.
	• MATLAB 6.0 or newer running.


3. Installation

	Unpack the archive in a folder that is part of the MATLAB path. The class
	functions reside in the '@extended' folder. For Pentium III compatible
	computers or newer, the folder '@extended/SSE' contains optimized func-
	tions. Move them into the '@extended' folder for use, or delete them.


4. Comments

	The functions should always reside in a subfolder called '@extended' to make
	sure they are not called for any data except of type 'extended'. In general,
	they do not check the data type of input arguments.

   Constants

	The class provides access to extended precision constants. For a built-in
	function, MATLAB calls always the double precision function whenever no
	input argument is given. To get an extended precision constant, use:

		pi(extended)		pi
		eps(extended)		1.08E-19
		ones(extended)		1
		realmax(extended)	1.18E+4932
		realmin(extended)	3.37E-4932
		zeros(extended)		0

		rand(extended)		uniformly distributed random number
		randn(extended)		normally distributed random number

   Limitations

	The class provides all arithmetic functions, the matrix multiplication,
	and most elementary functions (see also the optimized double precision
	class). Currently, none of the inverse matrix functions (division, in-
	version, decompositions) has been implemented. Wrapper scripts are
	supplied, which call MATLAB built-in functions.

   Performance

	Due to the native execution by the processor, the performance is affected
	merely by memory throughput. An extended value occupies 10bytes instead of
	8bytes. Memory access may be faster if the data were aligned at 16byte
	address boundaries. Though, the current implementation minimizes memory
	consumption by packing consecutive values closely together.

   Numerical accuracy

	In general, a statement of "inverseFunction(function(value))" produces "value"
	with a relative error of less than 1E-16. The roundoff error of about 1.1E-19
	leads to relatively important deviations in exponentiation. Note also that in
	particular addition/subtraction in the argument of a logarithm are critical
	operations due to a relative amplification of the rounding error. The logarithm
	itself works accurate over the full complex plane R x iR. For the sake of per-
	formance, the inverse transcendental functions are currently not implemented in
	that way. See also the summary about complex functions (available online).


5. Copyright

	This class is published as freeware. The author reserves the right to modify
	any of the contained files.

	You are allowed to distribute the functions as long as you deliver for free the
	entire package.

		Path		Files

		/		Readme.txt
		@extended/	_colonobj.m	abs.dll		acos.dll
				acosh.dll	all.m		and.m
				angle.dll	any.m		asin.dll
				asinh.dll	atan.dll	atan2.dll
				atanh.dll	cat.m		ceil.dll
				char.dll	cis.dll		cis.m
				colon.dll	complex.dll	conj.dll
				cos.dll		cosh.dll	ctranspose.dll
				diag.m		disp.m		display.m
				double.dll	end.m		eps.dll
				eq.dll		exp.dll		extended.dll
				extended.m	eye.m		fix.dll
				floor.dll	fprintf.m	ge.dll
				gt.dll		horzcat.m	imag.dll
				int16.dll	int32.dll	int8.dll
				ipermute.m	isempty.m	isequal.m
				isfinite.dll	isinf.dll	isnan.dll
				isnumeric.m	isreal.m	ldivide.dll
				le.dll		length.m	log.dll
				log10.dll	log2.dll	lt.dll
				minus.dll	mldivide.m	mod.dll
				mpower.m	mrdivide.m	mtimes.dll
				ndims.m		ne.dll		not.m
				num2str.m	numel.m		ones.dll
				or.m		permute.m	pi.dll
				plus.dll	pow2.dll	power.m
				prod.dll	rand.dll	rand.m
				randn.dll	randn.m		rdivide.dll
				real.dll	realmax.m	realmin.m
				rem.dll		repmat.m	reshape.m
				round.dll	sign.dll	sin.dll
				single.dll	sinh.dll	size.dll
				sprintf.m	sqrt.dll	struct.m
				subsasgn.m	subsindex.dll	subsref.m
				sum.dll		tan.dll		tanh.dll
				times.dll	transpose.dll	tril.m
				triu.m		true.dll	uint16.dll
				uint32.dll	uint8.dll	uminus.dll
				uplus.dll	value.m		vertcat.m
				xor.m		zeros.dll
		@extended/sse/	mtimes.dll	true.dll


6. Warranty

	Any warranty is strictly refused and you cannot anticipate any financial or
	technical support in case of malfunction or damage.

	Feedback and comments are welcome. I will try to track reported problems and
	fix bugs.


7. History

   • March 2, 2005
	Initial relase.

   • March 11, 2005
	Service relase. Bug fix in "power.m".

   • May 2, 2005
	Service release. Function "numel.m" and version information added.

   • May 24, 2005
	Service release due to an error report by Carlos López. Function "num2str.m"
	added and "subsasgn.m" fixed for empty inputs.

   • July 24, 2005
	Bug fix in "prod.dll": callback to MATLAB's "sum" replaced by "prod".

   • May 8, 2008
	Bug fixes for "prod.dll", "mtimes.dll" and "sse/true.dll".
	Optimized routine for calculating the exponential according to
	Agner Fog, "Optimizing subroutines in assembly language," at
	Copenhagen University College of Engineering, www.agner.org.


8. Download

	Optimized MATLAB functions are available online at my webpage.
	Summaries are also published at MATLAB central:

			http://www.mathworks.com/matlabcentral/


9. Trademarks

	MATLAB is a registered trademark of The MathWorks, Inc. Pentium is a
	registered trademark of Intel Corporation. Other product or brand names
	are trademarks or registered trademarks of their respective holders.

		________________________________________________________

			      EPFL © 2005-2008, Lausanne
				Webmaster • 8 May 2008
