[segment .text align=16]
global	_fvdot
global	_fvdots
global	_fvdott


;void fvdot(double* or, double* oi, const double* sr, const double* si, const double* tr, const double* ti, int n)
;
_fvdot:
	push ebp
	push ebx
	push edi
	push esi
	fninit
	mov ebp,[esp+44]		; n
	mov edx,[esp+40]		; ti
	mov ecx,[esp+36]		; tr
	mov ebx,[esp+32]		; si
	mov eax,[esp+28]		; sr
	mov esi,[esp+24]		; oi
	mov edi,[esp+20]		; or
	test esi,esi
       jz near	.4
	test ebx,ebx
       jz near	.2
	test edx,edx
       jz near	.3
.1:	fld qword [eax]			; complex
	fld qword [ecx]
	fld st0
	fmul st0,st2			; ac
	fld qword [ebx]
	fmul st2,st0			; bc
	fld qword [edx]
	fmul st4,st0			; ad
	fmulp st1,st0			; bd
	faddp st1,st0			; ac + bd
	fxch st1,st0
	fsubp st2,st0			; ad - bc
	fld qword [eax+8]
	fld qword [ecx+8]
	fld st0
	fmul st0,st2
	add eax,byte 24
	fld qword [ebx+8]
	fmul st2,st0
	add ecx,byte 24
	fld qword [edx+8]
	fmul st4,st0
	add ebx,byte 24
	fmulp st1,st0
	add edx,byte 24
	faddp st1,st0			; ac + bd
	faddp st3,st0
	fsubp st1,st0			; ad - bc
	faddp st2,st0
	fld qword [eax-8]
	fld qword [ecx-8]
	fld st0
	fmul st0,st2
	fld qword [ebx-8]
	fmul st2,st0
	fld qword [edx-8]
	fmul st4,st0
	add edi,byte 8
	fmulp st1,st0
	add esi,byte 8
	faddp st1,st0
	sub ebp,byte 1
	faddp st3,st0
	fsubp st1,st0
	faddp st2,st0
	fstp qword [edi-8]
	fstp qword [esi-8]
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [eax]			; s real
	fld qword [ecx]
	fmul st0,st1
	fld qword [edx]
	fmulp st2,st0
	fld qword [eax+8]
	fld qword [ecx+8]
	fmul st0,st1
	add eax,byte 24
	faddp st2,st0
	fld qword [edx+8]
	fmulp st1,st0
	add ecx,byte 24
	faddp st2,st0
	add edx,byte 24
	fld qword [eax-8]
	fld qword [ecx-8]
	fmul st0,st1
	add edi,byte 8
	faddp st2,st0
	add esi,byte 8
	fld qword [edx-8]
	fmulp st1,st0
	sub ebp,byte 1
	faddp st2,st0
	fstp qword [edi-8]
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; t real
	fld qword [eax]
	fmul st0,st1
	fld qword [ebx]
	fmulp st2,st0
	fld qword [ecx+8]
	fld qword [eax+8]
	fmul st0,st1
	add ecx,byte 24
	faddp st2,st0
	fld qword [ebx+8]
	fmulp st1,st0
	add eax,byte 24
	faddp st2,st0
	add ebx,byte 24
	fld qword [ecx-8]
	fld qword [eax-8]
	fmul st0,st1
	add edi,byte 8
	faddp st2,st0
	fld qword [ebx-8]
	fmulp st1,st0
	add esi,byte 8
	faddp st2,st0
	sub ebp,byte 1
	fstp qword [edi-8]
	fchs
	fstp qword [esi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fmul qword [ecx]
	fld qword [eax+8]
	fmul qword [ecx+8]
	add eax,byte 24
	add ecx,byte 24
	faddp st1,st0
	fld qword [eax-8]
	fmul qword [ecx-8]
	add edi,byte 8
	faddp st1,st0
	sub ebp,byte 1
	fstp qword [edi-8]
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn


align	4


_fvdots:
	push ebp
	push ebx
	push edi
	push esi
	fninit
	mov ebp,[esp+44]		; n
	mov edx,[esp+40]		; ti
	mov ecx,[esp+36]		; tr
	mov ebx,[esp+32]		; si
	mov eax,[esp+28]		; sr
	mov esi,[esp+24]		; oi
	mov edi,[esp+20]		; or
	test esi,esi
       jz near	.4
	test ebx,ebx
       jz near	.2
	test edx,edx
       jz near	.3
.1:	fld qword [eax]			; complex
	fld qword [ecx]
	fld st0
	fmul st0,st2			; ac
	fld qword [ebx]
	fmul st2,st0			; bc
	fld qword [edx]
	fmul st4,st0			; ad
	fmulp st1,st0			; bd
	faddp st1,st0			; ac + bd
	fxch st1,st0
	fsubp st2,st0			; ad - bc
	fld qword [eax+8]
	fld qword [ecx+8]
	fld st0
	fmul st0,st2
	fld qword [ebx+8]
	fmul st2,st0
	add ecx,byte 24
	fld qword [edx+8]
	fmul st4,st0
	fmulp st1,st0
	add edx,byte 24
	faddp st1,st0			; ac + bd
	faddp st3,st0
	fsubp st1,st0			; ad - bc
	faddp st2,st0
	fld qword [eax+16]
	fld qword [ecx-8]
	fld st0
	fmul st0,st2
	fld qword [ebx+16]
	fmul st2,st0
	fld qword [edx-8]
	fmul st4,st0
	add edi,byte 8
	fmulp st1,st0
	add esi,byte 8
	faddp st1,st0
	sub ebp,byte 1
	faddp st3,st0
	fsubp st1,st0
	faddp st2,st0
	fstp qword [edi-8]
	fstp qword [esi-8]
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [eax]			; s real
	fld qword [ecx]
	fmul st0,st1
	fld qword [edx]
	fmulp st2,st0
	fld qword [eax+8]
	fld qword [ecx+8]
	fmul st0,st1
	add ecx,byte 24
	faddp st2,st0
	fld qword [edx+8]
	fmulp st1,st0
	add edx,byte 24
	faddp st2,st0
	fld qword [eax+16]
	fld qword [ecx-8]
	fmul st0,st1
	add edi,byte 8
	faddp st2,st0
	add esi,byte 8
	fld qword [edx-8]
	fmulp st1,st0
	sub ebp,byte 1
	faddp st2,st0
	fstp qword [edi-8]
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; t real
	fld qword [eax]
	fmul st0,st1
	fld qword [ebx]
	fmulp st2,st0
	fld qword [ecx+8]
	fld qword [eax+8]
	fmul st0,st1
	add ecx,byte 24
	faddp st2,st0
	fld qword [ebx+8]
	fmulp st1,st0
	faddp st2,st0
	fld qword [ecx-8]
	fld qword [eax+16]
	fmul st0,st1
	add edi,byte 8
	faddp st2,st0
	fld qword [ebx+16]
	fmulp st1,st0
	add esi,byte 8
	faddp st2,st0
	sub ebp,byte 1
	fstp qword [edi-8]
	fchs
	fstp qword [esi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fmul qword [ecx]
	fld qword [eax+8]
	fmul qword [ecx+8]
	add ecx,byte 24
	faddp st1,st0
	fld qword [eax+16]
	fmul qword [ecx-8]
	add edi,byte 8
	faddp st1,st0
	sub ebp,byte 1
	fstp qword [edi-8]
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn


align	4

_fvdott:
	push ebp
	push ebx
	push edi
	push esi
	fninit
	mov ebp,[esp+44]		; n
	mov edx,[esp+40]		; ti
	mov ecx,[esp+36]		; tr
	mov ebx,[esp+32]		; si
	mov eax,[esp+28]		; sr
	mov esi,[esp+24]		; oi
	mov edi,[esp+20]		; or
	test esi,esi
       jz near	.4
	test ebx,ebx
       jz near	.2
	test edx,edx
       jz near	.3
.1:	fld qword [eax]			; complex
	fld qword [ecx]
	fld st0
	fmul st0,st2			; ac
	fld qword [ebx]
	fmul st2,st0			; bc
	fld qword [edx]
	fmul st4,st0			; ad
	fmulp st1,st0			; bd
	faddp st1,st0			; ac + bd
	fxch st1,st0
	fsubp st2,st0			; ad - bc
	fld qword [eax+8]
	fld qword [ecx+8]
	fld st0
	fmul st0,st2
	add eax,byte 24
	fld qword [ebx+8]
	fmul st2,st0
	fld qword [edx+8]
	fmul st4,st0
	add ebx,byte 24
	fmulp st1,st0
	faddp st1,st0			; ac + bd
	faddp st3,st0
	fsubp st1,st0			; ad - bc
	faddp st2,st0
	fld qword [eax-8]
	fld qword [ecx+16]
	fld st0
	fmul st0,st2
	fld qword [ebx-8]
	fmul st2,st0
	fld qword [edx+16]
	fmul st4,st0
	add edi,byte 8
	fmulp st1,st0
	add esi,byte 8
	faddp st1,st0
	sub ebp,byte 1
	faddp st3,st0
	fsubp st1,st0
	faddp st2,st0
	fstp qword [edi-8]
	fstp qword [esi-8]
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [eax]			; s real
	fld qword [ecx]
	fmul st0,st1
	fld qword [edx]
	fmulp st2,st0
	fld qword [eax+8]
	fld qword [ecx+8]
	fmul st0,st1
	add eax,byte 24
	faddp st2,st0
	fld qword [edx+8]
	fmulp st1,st0
	faddp st2,st0
	fld qword [eax-8]
	fld qword [ecx+16]
	fmul st0,st1
	add edi,byte 8
	faddp st2,st0
	add esi,byte 8
	fld qword [edx+16]
	fmulp st1,st0
	sub ebp,byte 1
	faddp st2,st0
	fstp qword [edi-8]
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; t real
	fld qword [eax]
	fmul st0,st1
	fld qword [ebx]
	fmulp st2,st0
	fld qword [ecx+8]
	fld qword [eax+8]
	fmul st0,st1
	faddp st2,st0
	fld qword [ebx+8]
	fmulp st1,st0
	add eax,byte 24
	faddp st2,st0
	add ebx,byte 24
	fld qword [ecx+16]
	fld qword [eax-8]
	fmul st0,st1
	add edi,byte 8
	faddp st2,st0
	fld qword [ebx-8]
	fmulp st1,st0
	add esi,byte 8
	faddp st2,st0
	sub ebp,byte 1
	fstp qword [edi-8]
	fchs
	fstp qword [esi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fmul qword [ecx]
	fld qword [eax+8]
	fmul qword [ecx+8]
	add eax,byte 24
	faddp st1,st0
	fld qword [eax-8]
	fmul qword [ecx+16]
	add edi,byte 8
	faddp st1,st0
	sub ebp,byte 1
	fstp qword [edi-8]
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
