
				     MATLAB toolbox
		________________________________________________________

				    Vector functions
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

	Marcel Leutenegger		marcel.leutenegger@epfl.ch
	EPFL STI IOA LOB
	BM 4.143			Phone:	+41 21 693 77 19
	Station 17
	CH-1015 Lausanne



1. Introduction

	The MATLAB toolbox offers two particular functions working on vectors as the
	vector cross product and the dot product. These functions are implemented to
	work on a general matrix dimension. Unfortunately, this implies a long overhead
	and slows down the performance significantly.

	This archive offers some particular vector functions as the vector norm, the
	vector length, scaling and the cross respectively dot product. All functions
	work on the first dimension of [3xM] matrices.


2. Requirements

	• An Intel Pentium II compatible computer or newer.
	• MATLAB 6.0 or newer running.


3. Installation

	Unpack the archive in a folder that is part of the MATLAB path. The subfolder
	'@class' contains the extension libraries.

	If you unpack the source code, the subfolder 'source' contains all necessary C
	and assembler files. Each function splits in a C-MEX wrapper part to be trans-
	lated with MATLAB:

		mex [-inline] [-O] func.c ffunc.obj [fdim.obj]

	The underlying assembler functions are written in the Netwide Assembler syntax.
	It is freely available on the Internet. The command is:

		nasmw -f win32 ffunc.asm

	Each function comes with an equivalent MATLAB script that also provides the
	command line help.


4. Comments

	The libraries should always reside in a subfolder called '@class' to make sure
	they are not called for any data except of type 'class'. They do not check the
	data type of passed arguments.

	You may consult the C-MEX wrapper performance hints respectively the floating-
	point assembler hints.


5. Copyright

	These routines are published as freeware. The author reserves the right to modify
	any of the contained files.

	You are allowed to distribute this package as long as you deliver the entire,
	original package for free.

		Path		File		Description

		/		Readme.txt	This manual
		@class/		vabs.dll	Vector length
				vabs.m		Help file
				vadd.dll	Vector addition
				vadd.m		Help file
				vcross.dll	Vector cross product
				vcross.m	Help file
				vdot.dll	Vector dot product
				vdot.m		Help file
				vmul.dll	Vector scaling
				vmul.m		Help file
				vnorm.dll	Vector norm
				vnorm.m		Help file

	The source code is fully disclosed to serve as practical example and for
	personal use only. You are not allowed to distribute the source code but you
	may apply it for your personal issues.

		Path		File		Description

		/ 		CmexWrapper.txt	C-MEX wrapper performance hints
				FPUasmHints.txt	Floating-point assembler hints
				Readme.txt	This manual
		source/		fdim.asm	Compare dimensions
				fvabs.asm	Vector length
				fvadd.asm	Vector addition
				fvcross.asm	Vector cross product
				fvdot.asm	Vector dot product
				fvmul.asm	Vector scaling
				fvnorm.asm	Vector norm
				vabs.c		Vector length
				vabs.m		Help file
				vadd.c		Vector addition
				vadd.m		Help file
				vcross.c	Vector cross product
				vcross.m	Help file
				vdot.c		Vector dot product
				vdot.m		Help file
				vmul.c		Vector scaling
				vmul.m		Help file
				vnorm.c		Vector norm
				vnorm.m		Help file


6. Warranty

	Any warranty is strictly refused and you cannot anticipate any financial or
	technical support in case of malfunction or damage.

	Feedback and comments are welcome. I will try to track reported problems and
	fix bugs.


7. History

   • January 18, 2004
	Initial relase.

   • March 3, 2005
	Service relase. Bug fix of the scripts.

   • March 6, 2005
	Initial extended precision release.

   • May 2, 2005
	Service release. Version information included.

   • December 29, 2005
	Initial single precision release and "vadd" included.


8. Download

	Optimized MATLAB routines are available online at

		   http://ioalinux1.epfl.ch/~mleutene/MATLABToolbox/


	A summary is also put for download at MATLAB central:

			http://www.mathworks.com/matlabcentral/


9. Trademarks

	MATLAB is a registered trademark of The MathWorks, Inc. Pentium is a
	registered trademark of Intel Corporation. Other product or brand names
	are trademarks or registered trademarks of their respective holders.

		________________________________________________________

			    Site map • EPFL © 2005, Lausanne
			      Webmaster • 29 December 2005
