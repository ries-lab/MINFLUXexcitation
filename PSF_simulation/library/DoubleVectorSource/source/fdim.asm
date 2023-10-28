[segment .text]
global	_fdim


;bool fdim(int m, int n, const int* d, const int* e)
;
_fdim:
	mov ecx,[esp+4]
	cmp ecx,[esp+8]
       jnz	.0
	push edi
	push esi
	pushfd
	cld
	mov edi,[esp+24]
	mov esi,[esp+28]
	repz cmpsd
       jnz	.1
	popfd
	pop esi
	pop edi
	mov eax,1
       retn
.1:	popfd
	pop esi
	pop edi
.0:	xor eax,eax
       retn
