[segment .text align=16]
global	_fvmul
global	_fvmuls
global	_fvmult


;void fvmul(double* or, double* oi, const double* sr, const double* si, const double* tr, const double* ti, int n)
;
_fvmul:
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
       jz near	.3
	test edx,edx
       jz	.2
.1:	fld qword [ecx]			; o complex
	fld qword [edx]
	fld qword [eax]
	fld qword [ebx]
	fld st3
	fmul st0,st2
	add ecx,byte 8
	fld st3
	fmul st0,st2
	add edx,byte 8
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st3
	fxch st1,st0
	fmul st0,st2
	faddp st1,st0
	fstp qword [esi]
	fld qword [eax+8]
	fld qword [ebx+8]
	fld st3
	fmul st0,st2
	add eax,byte 24
	fld st3
	fmul st0,st2
	add ebx,byte 24
	fsubp st1,st0
	fstp qword [edi+8]
	fmul st0,st3
	add edi,byte 24
	fxch st1,st0
	fmul st0,st2
	faddp st1,st0
	fstp qword [esi+8]
	fld qword [eax-8]
	fld qword [ebx-8]
	fld st3
	fmul st0,st2
	add esi,byte 24
	fld st3
	fmul st0,st2
	sub ebp,byte 1
	fsubp st1,st0
	fstp qword [edi-8]
	fmulp st3,st0
	fmulp st1,st0
	faddp st1,st0
	fstp qword [esi-8]
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [ecx]			; t real
	fld st0
	fmul qword [eax]
	add ecx,byte 8
	fstp qword [edi]
	fld st0
	fmul qword [ebx]
	fstp qword [esi]
	fld st0
	fmul qword [eax+8]
	add eax,byte 24
	fstp qword [edi+8]
	add edi,byte 24
	fld st0
	fmul qword [ebx+8]
	add ebx,byte 24
	fstp qword [esi+8]
	add esi,byte 24
	fld st0
	fmul qword [eax-8]
	sub ebp,byte 1
	fstp qword [edi-8]
	fmul qword [ebx-8]
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; s real
	fld qword [edx]
	fld qword [eax]
	fld st2
	fmul st0,st1
	add ecx,byte 8
	fstp qword [edi]
	fmul st0,st1
	add edx,byte 8
	fstp qword [esi]
	fld qword [eax+8]
	fld st2
	fmul st0,st1
	add eax,byte 24
	fstp qword [edi+8]
	fmul st0,st1
	add edi,byte 24
	fstp qword [esi+8]
	add esi,byte 24
	fld qword [eax-8]
	fmul st2,st0
	sub ebp,byte 1
	fmulp st1,st0
	fstp qword [esi-8]
	fstp qword [edi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [ecx]			; o real
	fld st0
	fmul qword [eax]
	add ecx,byte 8
	fstp qword [edi]
	fld st0
	fmul qword [eax+8]
	add eax,byte 24
	fstp qword [edi+8]
	fmul qword [eax-8]
	add edi,byte 24
	sub ebp,byte 1
	fstp qword [edi-8]
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn


align	4

_fvmuls:
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
       jz near	.3
	test edx,edx
       jz	.2
.1:	fld qword [ecx]			; complex
	fld qword [edx]
	fld qword [eax]
	fld qword [ebx]
	fld st3
	fmul st0,st2
	add ecx,byte 8
	fld st3
	fmul st0,st2
	add edx,byte 8
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st3
	fxch st1,st0
	fmul st0,st2
	faddp st1,st0
	fstp qword [esi]
	fld qword [eax+8]
	fld qword [ebx+8]
	fld st3
	fmul st0,st2
	fld st3
	fmul st0,st2
	fsubp st1,st0
	fstp qword [edi+8]
	fmul st0,st3
	add edi,byte 24
	fxch st1,st0
	fmul st0,st2
	faddp st1,st0
	fstp qword [esi+8]
	fld qword [eax+16]
	fld qword [ebx+16]
	fld st3
	fmul st0,st2
	add esi,byte 24
	fld st3
	fmul st0,st2
	sub ebp,byte 1
	fsubp st1,st0
	fstp qword [edi-8]
	fmulp st3,st0
	fmulp st1,st0
	faddp st1,st0
	fstp qword [esi-8]
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [ecx]			; t real
	fld qword [eax]
	fmul st0,st1
	add ecx,byte 8
	fstp qword [edi]
	fld qword [ebx]
	fmul st0,st1
	fstp qword [esi]
	fld qword [eax+8]
	fmul st0,st1
	fstp qword [edi+8]
	fld qword [ebx+8]
	fmul st0,st1
	add edi,byte 24
	fstp qword [esi+8]
	add esi,byte 24
	fld qword [eax+16]
	fmul st0,st1
	sub ebp,byte 1
	fstp qword [edi-8]
	fmul qword [ebx+16]
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; s real
	fld qword [edx]
	fld qword [eax]
	fld st0
	fmul st0,st3
	add ecx,byte 8
	fstp qword [edi]
	fmul st0,st1
	add edx,byte 8
	fstp qword [esi]
	fld qword [eax+8]
	fld st0
	fmul st0,st3
	fstp qword [edi+8]
	fmul st0,st1
	add edi,byte 24
	fstp qword [esi+8]
	add esi,byte 24
	fld qword [eax+16]
	fmul st2,st0
	sub ebp,byte 1
	fmulp st1,st0
	fstp qword [esi-8]
	fstp qword [edi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fmul qword [ecx]
	fstp qword [edi]
	fld qword [eax+8]
	fmul qword [ecx]
	fstp qword [edi+8]
	add edi,byte 24
	fld qword [eax+16]
	fmul qword [ecx]
	add ecx,byte 8
	sub ebp,byte 1
	fstp qword [edi-8]
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn


align	4

_fvmult:
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
       jz near	.3
	test edx,edx
       jz	.2
.1:	fld qword [ecx]			; complex
	fld qword [edx]
	fld qword [eax]
	fld qword [ebx]
	fld st3
	fmul st0,st2
	fld st3
	fmul st0,st2
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st3
	fxch st1,st0
	fmul st0,st2
	faddp st1,st0
	fstp qword [esi]
	fld qword [eax+8]
	fld qword [ebx+8]
	fld st3
	fmul st0,st2
	add eax,byte 24
	fld st3
	fmul st0,st2
	add ebx,byte 24
	fsubp st1,st0
	fstp qword [edi+8]
	fmul st0,st3
	add edi,byte 24
	fxch st1,st0
	fmul st0,st2
	faddp st1,st0
	fstp qword [esi+8]
	fld qword [eax-8]
	fld qword [ebx-8]
	fld st3
	fmul st0,st2
	add esi,byte 24
	fld st3
	fmul st0,st2
	sub ebp,byte 1
	fsubp st1,st0
	fstp qword [edi-8]
	fmulp st3,st0
	fmulp st1,st0
	faddp st1,st0
	fstp qword [esi-8]
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [ecx]			; t real
	fld st0
	fmul qword [eax]
	fstp qword [edi]
	fld st0
	fmul qword [eax+8]
	add eax,byte 24
	fstp qword [edi+8]
	fld st0
	fmul qword [ebx]
	add edi,byte 24
	fstp qword [esi]
	fld st0
	fmul qword [ebx+8]
	add ebx,byte 24
	fstp qword [esi+8]
	add esi,byte 24
	fld st0
	fmul qword [eax-8]
	sub ebp,byte 1
	fstp qword [edi-8]
	fmul qword [ebx-8]
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; s real
	fld qword [edx]
	fld qword [eax]
	fld st0
	fmul st0,st3
	fstp qword [edi]
	fmul st0,st1
	fstp qword [esi]
	fld qword [eax+8]
	fld st0
	fmul st0,st3
	add eax,byte 24
	fstp qword [edi+8]
	fmul st0,st1
	add edi,byte 24
	fstp qword [esi+8]
	add esi,byte 24
	fld qword [eax-8]
	fmul st2,st0
	sub ebp,byte 1
	fmulp st1,st0
	fstp qword [esi-8]
	fstp qword [edi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [ecx]			; o real
	fld st0
	fmul qword [eax]
	fstp qword [edi]
	fld st0
	fmul qword [eax+8]
	add eax,byte 24
	fstp qword [edi+8]
	add edi,byte 24
	fmul qword [eax-8]
	sub ebp,byte 1
	fstp qword [edi-8]
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
