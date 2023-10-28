[segment .text align=16]
global	_fvcross
global	_fvcrosss
global	_fvcrosst


;void fvcross(double* or, double* oi, const double* sr, const double* si, const double* tr, const double* ti, int n)
;
_fvcross:
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
	fld qword [eax+8]
	fld qword [ebx]
	fld qword [ebx+8]
	fld qword [eax+16]
	fld qword [ebx+16]
	fld st4
	fmul qword [ecx+16]
	add eax,byte 24
	fld st3
	fmul qword [edx+16]
	add ebx,byte 24
	fsubp st1,st0
	fld st2
	fmul qword [ecx+8]
	fsubp st1,st0
	fld st1
	fmul qword [edx+8]
	faddp st1,st0
	fstp qword [edi]		; Re(ox) = Re(sy)*Re(tz) - Im(sy)*Im(tz) - Re(sz)*Re(ty) + Im(sz)*Im(ty)
	fld st4
	fmul qword [edx+16]
	add edi,byte 24
	fld st3
	fmul qword [ecx+16]
	faddp st1,st0
	fld st2
	fmul qword [edx+8]
	fsubp st1,st0
	fld st1
	fmul qword [ecx+8]
	fsubp st1,st0
	fstp qword [esi]		; Im(ox) = Re(sy)*Im(tz) + Im(sy)*Re(tz) - Re(sz)*Im(ty) - Im(sz)*Im(ty)
	fld st1
	fmul qword [ecx]
	add esi,byte 24
	fld st1
	fmul qword [edx]
	fsubp st1,st0
	fld st6
	fmul qword [ecx+16]
	fsubp st1,st0
	fld st4
	fmul qword [edx+16]
	faddp st1,st0
	fstp qword [edi-16]		; Re(oy)
	fmul qword [ecx]
	fxch st1,st0
	fmul qword [edx]
	faddp st1,st0
	fld st4
	fmul qword [edx+16]
	fsubp st1,st0
	fld st2
	fmul qword [ecx+16]
	fsubp st1,st0
	fstp qword [esi-16]		; Im(oy)
	fld st3
	fmul qword [ecx+8]
	fld st2
	fmul qword [edx+8]
	fsubp st1,st0
	fld st3
	fmul qword [ecx]
	fsubp st1,st0
	fld st1
	fmul qword [edx]
	faddp st1,st0
	fstp qword [edi-8]		; Re(oz)
	fmul qword [ecx]
	fxch st1,st0
	fmul qword [ecx+8]
	add ecx,byte 24
	fsubrp st1,st0
	fxch st1,st0
	fmul qword [edx]
	fsubp st1,st0
	fxch st1,st0
	fmul qword [edx+8]
	add edx,byte 24
	faddp st1,st0
	sub ebp,byte 1
	fstp qword [esi-8]		; Im(oz)
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [eax]			; s real
	fld qword [eax+8]
	fld qword [eax+16]
	fld qword [ecx]
	fld qword [ecx+8]
	fld qword [ecx+16]
	fld st4
	fmul st0,st1
	add eax,byte 24
	fld st4
	fmul st0,st3
	add ecx,byte 24
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st5
	add edi,byte 24
	fld st3
	fmul st0,st3
	fsubrp st1,st0
	fstp qword [edi-16]
	fmul st0,st4
	fxch st1,st0
	fmul st0,st3
	fsubp st1,st0
	fstp qword [edi-8]
	fld qword [edx]
	fld qword [edx+8]
	fld qword [edx+16]
	add edx,byte 24
	fld st4
	fmul st0,st1
	fld st4
	fmul st0,st3
	fsubp st1,st0
	fstp qword [esi]
	fmul st0,st5
	add esi,byte 24
	fxch st3,st0
	fmul st0,st2
	sub ebp,byte 1
	fsubrp st3,st0
	fmulp st4,st0
	fmulp st2,st0
	fstp qword [esi-16]
	fsubp st1,st0
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; t real
	fld qword [ecx+8]
	fld qword [ecx+16]
	fld qword [eax]
	fld qword [eax+8]
	fld qword [eax+16]
	fld st1
	fmul st0,st4
	add ecx,byte 24
	fld st1
	fmul st0,st6
	add eax,byte 24
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st5
	add edi,byte 24
	fld st2
	fmul st0,st4
	fsubp st1,st0
	fstp qword [edi-16]
	fmul st0,st4
	fxch st1,st0
	fmul st0,st3
	fsubrp st1,st0
	fstp qword [edi-8]
	fld qword [ebx]
	fld qword [ebx+8]
	fld qword [ebx+16]
	add ebx,byte 24
	fld st1
	fmul st0,st4
	fld st1
	fmul st0,st6
	fsubp st1,st0
	fstp qword [esi]
	fmul st0,st5
	add esi,byte 24
	fxch st3,st0
	fmul st0,st2
	sub ebp,byte 1
	fsubp st3,st0
	fmulp st4,st0
	fmulp st2,st0
	fstp qword [esi-16]
	fsubrp st1,st0
	fstp qword [esi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fld qword [eax+8]
	fld qword [ecx]
	fld qword [ecx+8]
	fld qword [eax+16]
	fld qword [ecx+16]
	fld st4
	fmul st0,st1
	add eax,byte 24
	fld st2
	fmul st0,st4
	add ecx,byte 24
	fsubp st1,st0
	fstp qword [edi]		; ox
	fld st3
	fmulp st2,st0
	add edi,byte 24
	fmul st0,st5
	sub ebp,byte 1
	fsubp st1,st0
	fstp qword [edi-16]		; oy
	fmulp st3,st0
	fmulp st1,st0
	fsubp st1,st0
	fstp qword [edi-8]		; oz
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn


align	4

_fvcrosss:
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
	fld qword [eax+8]
	fld qword [ebx]
	fld qword [ebx+8]
	fld qword [eax+16]
	fld qword [ebx+16]
	fld st4
	fmul qword [ecx+16]
	fld st3
	fmul qword [edx+16]
	fsubp st1,st0
	fld st2
	fmul qword [ecx+8]
	fsubp st1,st0
	fld st1
	fmul qword [edx+8]
	faddp st1,st0
	fstp qword [edi]		; Re(ox) = Re(sy)*Re(tz) - Im(sy)*Im(tz) - Re(sz)*Re(ty) + Im(sz)*Im(ty)
	fld st4
	fmul qword [edx+16]
	add edi,byte 24
	fld st3
	fmul qword [ecx+16]
	faddp st1,st0
	fld st2
	fmul qword [edx+8]
	fsubp st1,st0
	fld st1
	fmul qword [ecx+8]
	fsubp st1,st0
	fstp qword [esi]		; Im(ox) = Re(sy)*Im(tz) + Im(sy)*Re(tz) - Re(sz)*Im(ty) - Im(sz)*Im(ty)
	fld st1
	fmul qword [ecx]
	add esi,byte 24
	fld st1
	fmul qword [edx]
	fsubp st1,st0
	fld st6
	fmul qword [ecx+16]
	fsubp st1,st0
	fld st4
	fmul qword [edx+16]
	faddp st1,st0
	fstp qword [edi-16]		; Re(oy)
	fmul qword [ecx]
	fxch st1,st0
	fmul qword [edx]
	faddp st1,st0
	fld st4
	fmul qword [edx+16]
	fsubp st1,st0
	fld st2
	fmul qword [ecx+16]
	fsubp st1,st0
	fstp qword [esi-16]		; Im(oy)
	fld st3
	fmul qword [ecx+8]
	fld st2
	fmul qword [edx+8]
	fsubp st1,st0
	fld st3
	fmul qword [ecx]
	fsubp st1,st0
	fld st1
	fmul qword [edx]
	faddp st1,st0
	fstp qword [edi-8]		; Re(oz)
	fmul qword [ecx]
	fxch st1,st0
	fmul qword [ecx+8]
	add ecx,byte 24
	fsubrp st1,st0
	fxch st1,st0
	fmul qword [edx]
	fsubp st1,st0
	fxch st1,st0
	fmul qword [edx+8]
	add edx,byte 24
	faddp st1,st0
	sub ebp,byte 1
	fstp qword [esi-8]		; Im(oz)
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [eax]			; s real
	fld qword [eax+8]
	fld qword [eax+16]
	fld qword [ecx]
	fld qword [ecx+8]
	fld qword [ecx+16]
	fld st4
	fmul st0,st1
	add ecx,byte 24
	fld st4
	fmul st0,st3
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st5
	add edi,byte 24
	fld st3
	fmul st0,st3
	fsubrp st1,st0
	fstp qword [edi-16]
	fmul st0,st4
	fxch st1,st0
	fmul st0,st3
	fsubp st1,st0
	fstp qword [edi-8]
	fld qword [edx]
	fld qword [edx+8]
	fld qword [edx+16]
	add edx,byte 24
	fld st4
	fmul st0,st1
	fld st4
	fmul st0,st3
	fsubp st1,st0
	fstp qword [esi]
	fmul st0,st5
	add esi,byte 24
	fxch st3,st0
	fmul st0,st2
	sub ebp,byte 1
	fsubrp st3,st0
	fmulp st4,st0
	fmulp st2,st0
	fstp qword [esi-16]
	fsubp st1,st0
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; t real
	fld qword [ecx+8]
	fld qword [ecx+16]
	fld qword [eax]
	fld qword [eax+8]
	fld qword [eax+16]
	fld st1
	fmul st0,st4
	add ecx,byte 24
	fld st1
	fmul st0,st6
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st5
	add edi,byte 24
	fld st2
	fmul st0,st4
	fsubp st1,st0
	fstp qword [edi-16]
	fmul st0,st4
	fxch st1,st0
	fmul st0,st3
	fsubrp st1,st0
	fstp qword [edi-8]
	fld qword [ebx]
	fld qword [ebx+8]
	fld qword [ebx+16]
	fld st1
	fmul st0,st4
	fld st1
	fmul st0,st6
	fsubp st1,st0
	fstp qword [esi]
	fmul st0,st5
	add esi,byte 24
	fxch st3,st0
	fmul st0,st2
	sub ebp,byte 1
	fsubp st3,st0
	fmulp st4,st0
	fmulp st2,st0
	fstp qword [esi-16]
	fsubrp st1,st0
	fstp qword [esi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fld qword [eax+8]
	fld qword [ecx]
	fld qword [ecx+8]
	fld qword [eax+16]
	fld qword [ecx+16]
	fld st4
	fmul st0,st1
	add ecx,byte 24
	fld st2
	fmul st0,st4
	fsubp st1,st0
	fstp qword [edi]		; ox
	fld st3
	fmulp st2,st0
	add edi,byte 24
	fmul st0,st5
	sub ebp,byte 1
	fsubp st1,st0
	fstp qword [edi-16]		; oy
	fmulp st3,st0
	fmulp st1,st0
	fsubp st1,st0
	fstp qword [edi-8]		; oz
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn


align	4

_fvcrosst:
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
	fld qword [eax+8]
	fld qword [ebx]
	fld qword [ebx+8]
	fld qword [eax+16]
	fld qword [ebx+16]
	fld st4
	fmul qword [ecx+16]
	add eax,byte 24
	fld st3
	fmul qword [edx+16]
	add ebx,byte 24
	fsubp st1,st0
	fld st2
	fmul qword [ecx+8]
	fsubp st1,st0
	fld st1
	fmul qword [edx+8]
	faddp st1,st0
	fstp qword [edi]		; Re(ox) = Re(sy)*Re(tz) - Im(sy)*Im(tz) - Re(sz)*Re(ty) + Im(sz)*Im(ty)
	fld st4
	fmul qword [edx+16]
	add edi,byte 24
	fld st3
	fmul qword [ecx+16]
	faddp st1,st0
	fld st2
	fmul qword [edx+8]
	fsubp st1,st0
	fld st1
	fmul qword [ecx+8]
	fsubp st1,st0
	fstp qword [esi]		; Im(ox) = Re(sy)*Im(tz) + Im(sy)*Re(tz) - Re(sz)*Im(ty) - Im(sz)*Im(ty)
	fld st1
	fmul qword [ecx]
	add esi,byte 24
	fld st1
	fmul qword [edx]
	fsubp st1,st0
	fld st6
	fmul qword [ecx+16]
	fsubp st1,st0
	fld st4
	fmul qword [edx+16]
	faddp st1,st0
	fstp qword [edi-16]		; Re(oy)
	fmul qword [ecx]
	fxch st1,st0
	fmul qword [edx]
	faddp st1,st0
	fld st4
	fmul qword [edx+16]
	fsubp st1,st0
	fld st2
	fmul qword [ecx+16]
	fsubp st1,st0
	fstp qword [esi-16]		; Im(oy)
	fld st3
	fmul qword [ecx+8]
	fld st2
	fmul qword [edx+8]
	fsubp st1,st0
	fld st3
	fmul qword [ecx]
	fsubp st1,st0
	fld st1
	fmul qword [edx]
	faddp st1,st0
	fstp qword [edi-8]		; Re(oz)
	fmul qword [ecx]
	sub ebp,byte 1
	fxch st1,st0
	fmul qword [ecx+8]
	fsubrp st1,st0
	fxch st1,st0
	fmul qword [edx]
	fsubp st1,st0
	fxch st1,st0
	fmul qword [edx+8]
	faddp st1,st0
	fstp qword [esi-8]		; Im(oz)
       jg	.1
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.2:	fld qword [eax]			; s real
	fld qword [eax+8]
	fld qword [eax+16]
	fld qword [ecx]
	fld qword [ecx+8]
	fld qword [ecx+16]
	fld st4
	fmul st0,st1
	add eax,byte 24
	fld st4
	fmul st0,st3
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st5
	add edi,byte 24
	fld st3
	fmul st0,st3
	fsubrp st1,st0
	fstp qword [edi-16]
	fmul st0,st4
	fxch st1,st0
	fmul st0,st3
	fsubp st1,st0
	fstp qword [edi-8]
	fld qword [edx]
	fld qword [edx+8]
	fld qword [edx+16]
	fld st4
	fmul st0,st1
	fld st4
	fmul st0,st3
	fsubp st1,st0
	fstp qword [esi]
	fmul st0,st5
	add esi,byte 24
	fxch st3,st0
	fmul st0,st2
	sub ebp,byte 1
	fsubrp st3,st0
	fmulp st4,st0
	fmulp st2,st0
	fstp qword [esi-16]
	fsubp st1,st0
	fstp qword [esi-8]
       jg	.2
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.3:	fld qword [ecx]			; t real
	fld qword [ecx+8]
	fld qword [ecx+16]
	fld qword [eax]
	fld qword [eax+8]
	fld qword [eax+16]
	fld st1
	fmul st0,st4
	add eax,byte 24
	fld st1
	fmul st0,st6
	fsubp st1,st0
	fstp qword [edi]
	fmul st0,st5
	add edi,byte 24
	fld st2
	fmul st0,st4
	fsubp st1,st0
	fstp qword [edi-16]
	fmul st0,st4
	fxch st1,st0
	fmul st0,st3
	fsubrp st1,st0
	fstp qword [edi-8]
	fld qword [ebx]
	fld qword [ebx+8]
	fld qword [ebx+16]
	fld st1
	fmul st0,st4
	add ebx,byte 24
	fld st1
	fmul st0,st6
	fsubp st1,st0
	fstp qword [esi]
	fmul st0,st5
	add esi,byte 24
	fxch st3,st0
	fmul st0,st2
	sub ebp,byte 1
	fsubp st3,st0
	fmulp st4,st0
	fmulp st2,st0
	fstp qword [esi-16]
	fsubrp st1,st0
	fstp qword [esi-8]
       jg	.3
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
.4:	fld qword [eax]			; o real
	fld qword [eax+8]
	fld qword [ecx]
	fld qword [ecx+8]
	fld qword [eax+16]
	fld qword [ecx+16]
	fld st4
	fmul st0,st1
	add eax,byte 24
	fld st2
	fmul st0,st4
	fsubp st1,st0
	fstp qword [edi]		; ox
	fld st3
	fmulp st2,st0
	add edi,byte 24
	fmul st0,st5
	sub ebp,byte 1
	fsubp st1,st0
	fstp qword [edi-16]		; oy
	fmulp st3,st0
	fmulp st1,st0
	fsubp st1,st0
	fstp qword [edi-8]		; oz
       jg	.4
	pop esi
	pop edi
	pop ebx
	pop ebp
       retn
