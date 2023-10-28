[segment .text align=16]
global	_fvabs


;void fvabs(double* or, const double* sr, const double* si, int n)
;
_fvabs:	push ebx
	fninit
	mov ecx,[esp+20]	; n
	mov ebx,[esp+16]	; si
	mov eax,[esp+12]	; sr
	test ebx,ebx
	mov edx,[esp+8]		; or
       jz	.2
.1:	fld qword [eax]
	fmul st0,st0
	fld qword [ebx]
	fmul st0,st0
	add edx,byte 8
	faddp st1,st0
	fld qword [eax+8]
	fmul st0,st0
	faddp st1,st0
	fld qword [ebx+8]
	fmul st0,st0
	faddp st1,st0
	fld qword [eax+16]
	fmul st0,st0
	add eax,byte 24
	faddp st1,st0
	fld qword [ebx+16]
	fmul st0,st0
	add ebx,byte 24
	faddp st1,st0
	sub ecx,byte 1
	fsqrt
	fstp qword [edx-8]
       jg	.1
	pop ebx
       retn
.2:	fld qword [eax]
	fmul st0,st0
	fld qword [eax+8]
	fmul st0,st0
	add edx,byte 8
	faddp st1,st0
	fld qword [eax+16]
	fmul st0,st0
	add eax,byte 24
	faddp st1,st0
	sub ecx,byte 1
	fsqrt
	fstp qword [edx-8]
       jg	.2
	pop ebx
       retn
