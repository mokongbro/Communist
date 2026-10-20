	.file	"trans.cpp"
	.text
	.section	.text$_Z6printfPKcz,"x"
	.linkonce discard
	.p2align 4
	.globl	_Z6printfPKcz
	.def	_Z6printfPKcz;	.scl	2;	.type	32;	.endef
	.seh_proc	_Z6printfPKcz
_Z6printfPKcz:
.LFB9:
	pushq	%r12
	.seh_pushreg	%r12
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$56, %rsp
	.seh_stackalloc	56
	.seh_endprologue
	leaq	88(%rsp), %rbx
	movq	%rcx, %r12
	movq	%rdx, 88(%rsp)
	movl	$1, %ecx
	movq	%r8, 96(%rsp)
	movq	%r9, 104(%rsp)
	movq	%rbx, 40(%rsp)
	call	*__imp___acrt_iob_func(%rip)
	movq	%rbx, %r8
	movq	%r12, %rdx
	movq	%rax, %rcx
	call	__mingw_vfprintf
	addq	$56, %rsp
	popq	%rbx
	popq	%r12
	ret
	.seh_endproc
	.def	__main;	.scl	2;	.type	32;	.endef
	.section .rdata,"dr"
.LC0:
	.ascii "%s\12\0"
.LC1:
	.ascii "r\0"
.LC2:
	.ascii "w\0"
	.section	.text.startup,"x"
	.p2align 4
	.globl	main
	.def	main;	.scl	2;	.type	32;	.endef
	.seh_proc	main
main:
.LFB81:
	pushq	%rdi
	.seh_pushreg	%rdi
	pushq	%rsi
	.seh_pushreg	%rsi
	pushq	%rbx
	.seh_pushreg	%rbx
	subq	$32, %rsp
	.seh_stackalloc	32
	.seh_endprologue
	movl	%ecx, %esi
	movq	%rdx, %rbx
	call	__main
	cmpl	$5, %esi
	je	.L4
	testl	%esi, %esi
	jle	.L5
	leal	-1(%rsi), %eax
	leaq	8(%rbx,%rax,8), %rdi
	leaq	.LC0(%rip), %rsi
	.p2align 4,,10
	.p2align 3
.L6:
	movq	(%rbx), %rdx
	movq	%rsi, %rcx
	addq	$8, %rbx
	call	_Z6printfPKcz
	cmpq	%rdi, %rbx
	jne	.L6
.L5:
	xorl	%eax, %eax
	addq	$32, %rsp
	popq	%rbx
	popq	%rsi
	popq	%rdi
	ret
.L4:
	movq	__imp___acrt_iob_func(%rip), %rsi
	xorl	%ecx, %ecx
	call	*%rsi
	movq	8(%rbx), %rcx
	leaq	.LC1(%rip), %rdx
	movq	%rax, %r8
	call	freopen
	movl	$1, %ecx
	call	*%rsi
	movq	16(%rbx), %rcx
	leaq	.LC2(%rip), %rdx
	movq	%rax, %r8
	call	freopen
	movq	24(%rbx), %rcx
	call	atoi
	movq	32(%rbx), %rcx
	movl	%eax, %esi
	call	atoi
	movl	%eax, %ebx
	call	getchar
	cmpl	$-1, %eax
	movl	%eax, %ecx
	je	.L5
	.p2align 4,,10
	.p2align 3
.L9:
	cmpl	%ecx, %esi
	jne	.L17
	movl	%ebx, %ecx
.L17:
	call	putchar
	call	getchar
	cmpl	$-1, %eax
	movl	%eax, %ecx
	jne	.L9
	jmp	.L5
	.seh_endproc
	.ident	"GCC: (x86_64-posix-seh, Built by MinGW-Builds project) 11.4.0"
	.def	__mingw_vfprintf;	.scl	2;	.type	32;	.endef
	.def	freopen;	.scl	2;	.type	32;	.endef
	.def	atoi;	.scl	2;	.type	32;	.endef
	.def	getchar;	.scl	2;	.type	32;	.endef
	.def	putchar;	.scl	2;	.type	32;	.endef
