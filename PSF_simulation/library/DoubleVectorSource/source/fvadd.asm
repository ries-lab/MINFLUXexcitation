[segment .text align=16]
global	_fvadd
global	_fvadds
global	_fvmov
global	_fvmovs


;void fvadd(double* o, const double* s, const double* t, int n)
;
_fvadd:
	push ebx
	fninit
	mov ecx,[esp+20]		; n
	mov ebx,[esp+16]		; t
	mov eax,[esp+12]		; s
	mov edx,[esp+8]			; o
.1:	fld qword [eax]
	fadd qword [ebx]
	fld qword [eax+8]
	fadd qword [ebx+8]
	fld qword [eax+16]
	add eax,byte 24
	fadd qword [ebx+16]
	add ebx,byte 24
	fstp qword [edx+16]
	sub ecx,byte 1
	fstp qword [edx+8]
	fstp qword [edx]
	lea edx,[edx+24]
       jg	.1
	pop ebx
       retn


align	4

_fvadds:
	fninit
	mov edx,[esp+8]			; s
	mov eax,[esp+12]		; t
	fld qword [edx]
	mov ecx,[esp+16]		; n
	fld qword [edx+8]
	fld qword [edx+16]
	mov edx,[esp+4]			; o
.1:	fld qword [eax]
	fadd st0,st3
	fld qword [eax+8]
	fadd st0,st3
	fld qword [eax+16]
	fadd st0,st3
	add eax,byte 24
	fstp qword [edx+16]
	sub ecx,byte 1
	fstp qword [edx+8]
	fstp qword [edx]
	lea edx,[edx+24]
       jg	.1
	fninit
       retn


align	4

;void fvmov(double* o, const double* s, int n)
;
_fvmov:
	fninit
	mov ecx,[esp+12]		; n
	mov eax,[esp+8]			; s
	mov edx,[esp+4]			; o
.1:	fld qword [eax]
	fld qword [eax+8]
	fld qword [eax+16]
	add eax,byte 24
	fstp qword [edx+16]
	sub ecx,byte 1
	fstp qword [edx+8]
	fstp qword [edx]
	lea edx,[edx+24]
       jg	.1
       retn


align	4

_fvmovs:
	fninit
	mov ecx,[esp+12]		; n
	mov eax,[esp+8]			; s
	mov edx,[esp+4]			; o
	fld qword [eax+16]
	sub ecx,byte 1
	fld qword [eax+8]
	fld qword [eax]
       jng	.2
.1:	fst qword [edx]
	fld st1
	sub ecx,byte 1
	fstp qword [edx+8]
	fld st2
	fstp qword [edx+16]
	lea edx,[edx+24]
       jg	.1
.2:	fstp qword [edx]
	fstp qword [edx+8]
	fstp qword [edx+16]
       retn
