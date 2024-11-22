### A Pluto.jl notebook ###
# v0.20.0

#> [frontmatter]
#> title = "TETM241 | MATH557 | NOTES"
#> 
#>     [[frontmatter.author]]
#>     name = "Dr. Mohammed Alshahrani"
#>     url = "https://mshahrani.website/"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4eb18bb0-5b04-11ef-0c2c-8747a3f06685
begin
    using CommonMark
    using PlutoUI, PlutoExtras
    using Plots, PlotThemes, LaTeXStrings
    using Latexify
    using HypertextLiteral
    using Colors
    using LinearAlgebra, Random, Printf, SparseArrays
    using Symbolics, SymbolicUtils
    using QRCoders
    using PrettyTables
	using Combinatorics
	using Roots
	# using SymPy as sp
    # using NonlinearSolve
    # using ForwardDiff
    # using Integrals
	# using OrdinaryDiffEq
	# using DifferentialEquations
	# using ModelingToolkit
	
end

# ╔═╡ 01183bde-075c-4848-995a-b23ffeeb97c8
md"## mathmatize"

# ╔═╡ 73b56b54-22a2-4fa0-8eed-ab8a23cebc74
md"# Programming Assignments With Julia"

# ╔═╡ 06de9303-589d-40cf-ac95-4df2020af3a6
cm"""
#### Numerical Computation
In this course, implementation and programming assignments will be carried out through Julia programming language

##### Install Julia
Download and install Julia from [HERE](https://julialang.org/downloads/)

##### Learn Julia
1.	Ben Lauwens and Allen B. Downey. [Think Julia. O’Reilly Media, June 2019.](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html)
2.	[Other Books](https://julialang.org/learning/books/).


"""

# ╔═╡ c59dffb4-c419-461b-8096-e27171be0a87
# cm"""
# - LU and QR factorizations 2.4, 2.5, 3.2, 3.3, 3.5, 4.1, 4.2, 5.1− 5.4, 5.6
# - SVD, norms and LSQ 6.1, 6.3, 7.1 − 7.4, 8.1 − 8.3, 9.1 − 9.3, 9.4.1
# - Kronecker products 10.1, 10.2, 10.3, 11.1, 11.2, 11.3
# - Iterative methods 12.1 − 12.4, 13.1 − 13.3, 13.5
# - Eigenpairs 14.1 − 14.5, 15.1 − 15.3
# """

# ╔═╡ a8948e17-2846-431f-9765-6359eaeb20a9
md" # Chapter 1: Review of Basic Concepts."

# ╔═╡ db84f278-2b61-40fd-b0a8-bb132cff5f18
cm"""
> You should read this chapter to recall the basic concepts of linear algebra you learned in your undergarduate studies.
"""

# ╔═╡ df6c46d3-fa32-4709-a392-463167b33c46
md"## Column Space, Row Sapce, Null Space, Inner Product and Norm"

# ╔═╡ b976040e-4974-4bfc-a6a9-e3926a1f2eef
# using LinearAlgebra

# ╔═╡ 5356f40b-2cc7-4490-9752-115cae126839
let 
	X = [
	1 2 3 4 1
	1 6 6 9 4
	1 5 6 3 9]
	rank(X)
	# nullspace(X)
end

# ╔═╡ e3e2f379-f7bb-4ccf-85e0-378ff479f3de
1+2im

# ╔═╡ 6ccfccef-09cc-42d1-a5c8-3899801b4438


# ╔═╡ 17263bb2-6964-48a6-b31e-b20f168f35a2
# using LinearAlgebra

# ╔═╡ 30bc89ef-58af-4bed-b172-aa6b4e2c8491
let
	X = [1 1 2 1; 2 3 5 2; -2 2 1 -3]
	# x₁,x₂,x₃,x₄ = eachcol(X)
	# combs = combinations([1,2,3,4],3)
	# for a in combs
	# 	@show det(X[:,a])
	# end
	# nullspace(X)
	# Y = X'
	# nullspace(Y)
	# Z = [1+im 1 2 1;2 3 5 2;-2 2 1 -3]
	
end

# ╔═╡ 7e88c07e-394e-46eb-b1e5-89fb65287e36
md"##  Linear Systems"

# ╔═╡ 21a9d1b8-fdb4-44e6-99da-f97ac172a9a3
cm"""
Consider a linear system
```math
\begin{array}{cc}
a_{11} x_1+a_{12} x_2+\cdots+a_{1 n} x_n & =b_1 \\
a_{21} x_1+a_{22} x_2+\cdots+a_{2 n} x_n & =b_2 \\
\vdots & \vdots \\
\vdots & \vdots \\
a_{m 1} x_1+a_{m 2} x_2+\cdots+a_{m n} x_n & =b_m
\end{array}
```
This can be written
```math
\boldsymbol{A} \boldsymbol{x}=\left[\begin{array}{cccc}a_{11} & a_{12} & \cdots & a_{1 n} \\ a_{21} & a_{22} & \cdots & a_{2 n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m 1} & a_{m 2} & \cdots & a_{m n}\end{array}\right]\left[\begin{array}{c}x_1 \\ x_2 \\ \vdots \\ x_n\end{array}\right]=\left[\begin{array}{c}b_1 \\ b_2 \\ \vdots \\ b_m\end{array}\right]=\boldsymbol{b}
```

"""

# ╔═╡ c7f1e8cd-da40-469f-b8cf-42869d5b46ed
md"## The Inverse Matrix"

# ╔═╡ ffcc0c43-bbb5-433c-a98f-b8421a846485
md"## Determinants"

# ╔═╡ 8aa563a5-8264-4bf0-89c9-c1daa74bd4d6
cm"""
For any ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` the determinant of ``\boldsymbol{A}`` is defined by the number
```math
\operatorname{det}(\boldsymbol{A})=\sum_{\sigma \in S_n} \operatorname{sign}(\sigma) a_{\sigma(1), 1} a_{\sigma(2), 2} \cdots a_{\sigma(n), n}
```
"""

# ╔═╡ dbf3b4de-06ce-46ff-9178-cc28961ab3e5
function mydet(A)
	n,=size(A)
	σ= permutations(1:n,n)
	map(x->(x,parity(x)),σ) |> d -> map(x->(-1)^(x[2])*prod(A[x[1][i],i] for i in 1:n),d) |> sum
end

# ╔═╡ 5ab2ce16-f52d-4107-9035-4dc47df19fcd
let
	Random.seed!(0)
	n = 3
	A =[1 1 2 ; 2 3 5; -2 2 1] 
	# inv(A), mydet(A)
	# det(A)
	# # A = rand(-2:5,n,n)
	# σ= permutations(1:n,n)|>collect
	# [(s,prod(A[s[i],i] for i in 1:n)) for s in σ]
	# mydet = map(x->(x,parity(x)),σ) |> d -> map(x->(-1)^(x[2])*prod(A[x[1][i],i] for i in 1:n),d) |> sum
	# det(A),mydet
	# parity([3,5,1,2,4])
end

# ╔═╡ 4dbfa5bd-dfcc-4195-8337-02f8bed8748a
let
	Random.seed!(123)
	A= zeros(10,10)
	CD=A[1:5,:]=rand(2:8,5,10)
	C = CD[:,1:5]
	E=A[6:end,6:end]=rand(2:8,5,5)
	det(A),det(C[1:5,1:5])*det(E)
end

# ╔═╡ 49ee7daf-8b09-4303-8aa6-05ae57977688
let 
	# example
	A = [2 1;-1 3]
	# eigen(A)
	x=[1;-1]
	b = A*x
	A,b
end
	

# ╔═╡ 21e63aa1-c0de-4980-b973-9afd6b1dde87
md"## Eigenvalues, Eigenvectors and Eigenpairs"

# ╔═╡ 65875093-4ff0-4684-a355-678d727f36f7
cm"""
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is a square matrix, ``\lambda \in \mathbb{C}`` and ``\boldsymbol{x} \in \mathbb{C}^n``. 

- We say that ``(\lambda, \boldsymbol{x})`` is an eigenpair for ``\boldsymbol{A}`` if 
```math
\boldsymbol{A} \boldsymbol{x}=\lambda \boldsymbol{x} \quad \text{ and } \boldsymbol{x} \quad \text{is nonzero.}
``` 
- The scalar ``\lambda`` is called an __eigenvalue__ and ``\boldsymbol{x}`` is said to be an __eigenvector__. 
- The set of eigenvalues is called the __spectrum__ of ``\boldsymbol{A}`` and is denoted by ``\sigma(\boldsymbol{A})``. (For example, ``\sigma(\boldsymbol{I})=\{1, \ldots, 1\}=\{1\}``.)
-  For any ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` we have ``\lambda \in`` ``\sigma(\boldsymbol{A}) \Longleftrightarrow \operatorname{det}(\boldsymbol{A}-\lambda \boldsymbol{I})=0``
"""

# ╔═╡ 6d79d33f-8caa-4fa9-b1e3-baf6bb25d519
cm"""
- For any ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``
```math
\operatorname{trace}(\boldsymbol{A})=\lambda_1+\lambda_2+\cdots+\lambda_n, \quad \operatorname{det}(\boldsymbol{A})=\lambda_1 \lambda_2 \cdots \lambda_n
```
<div class="p40" >

where the trace of ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is the sum of its diagonal elements
```math
\operatorname{trace}(\boldsymbol{A}):=a_{11}+a_{22}+\cdots+a_{n n}
```
</div>

- The matrix ``A`` is singular if and only if zero is an eigenvalue.
"""

# ╔═╡ 374034ce-15d5-403d-9d2a-aa760cc0a269
md"# Chapter 2: Diagonally Dominant Tridiagonal Matrices"

# ╔═╡ e06e478f-b2a9-4b33-8a36-3bc87e53f62a
md"##  2.2 A Two Point Boundary Value Problem"

# ╔═╡ 8d801c77-08e3-45c0-b27c-0d36f0bfb29e
cm"""
Consider the simple two point boundary value problem

$(texeq"-u^{\prime \prime}(x)=f(x), \quad x \in[0,1], \quad u(0)=0, u(1)=0\label{1dpp}\tag{1-D Poisson}")

"""

# ╔═╡ 41e2ac3f-5488-4506-9278-b7857d8ad4e7
md"## LU Factorization of a Tridiagonal System"

# ╔═╡ d5e9beb5-56a9-4a99-903d-719dc18cf710
cm"""
Consider a tridiagonal linear system ``Ax = b`` where ``A = \operatorname{tridiag}(a_i,d_i,c_i) \in \mathbb{C}^{n\times n}``.

Construct ``L`` and ``U`` such that ``A=LU``. Then we have 
$(texeq"\left[\begin{array}{lllll}d_1 & c_1 & & & \\ a_1 & d_2 & c_2 & & \\ & \ddots & \ddots & \ddots & \\ & & a_{n-2} & d_{n-1} & c_{n-1} \\ & & & a_{n-1} & d_n\end{array}\right]=\left[\begin{array}{cccc}1 & & & \\ l_1 & 1 & \\ & \ddots & \ddots & \\ & & l_{n-1} & 1\end{array}\right]\left[\begin{array}{cccc}u_1 & c_1 & & \\ & \ddots & \ddots & \\ & & u_{n-1} & c_{n-1} \\ & & & u_n\end{array}\right]\label{lufactor}")
"""

# ╔═╡ 9a95f3c9-0e47-422a-9b0e-35fc8c528308
md"##  2.4 The Eigen pairs of the 1D Test Matrix"

# ╔═╡ 918ca6c0-ac2d-4072-b7d9-b73de95248fd
cm"""
The second derivative matrix ``\boldsymbol{T}=\operatorname{tridiag}(-1,2,-1)`` is a special case of the tridiagonal matrix
```math
\boldsymbol{T}_1:=\operatorname{tridiag}(a, d, a)
```
where ``a, d \in \mathbb{R}``. We call this the 1D test matrix. 
- It is symmetric and 
- strictly diagonally dominant if ``|d|>2|a|``.
- the eigenvectors are the columns of the sine matrix defined by
```math
\boldsymbol{S}=\left[\sin \frac{j k \pi}{m+1}\right]_{j, k=1}^m \in \mathbb{R}^{m \times m}
```
- the eigenvalues are given by
```math
\lambda_j=d+2 a \cos (j \pi /(m+1)), \quad j=1,\cdots, m.
```
- The eigenvalues of a Hermitian matrix are real. Moreover, eigenvectors corresponding to distinct eigen values are orthogonal.
"""

# ╔═╡ a5240ee4-1f4d-4229-95f4-fe115645c92c
begin
	triadiag(a,d,m) =diagm(0=>repeat([d],m),-1=>repeat([a],m-1),1=>repeat([a],m-1))
end

# ╔═╡ 05dfaf81-8b67-43de-b249-f1e23c3664f2
let
	a,d,m = -1,3,3
	A = triadiag(a,d,m)
	E,V = eigen(A) # using LinearAlgebra
	E2 = [d + 2*a*cos(j*π/(m+1)) for j in 1:m]
	V1= [sin(j*k*π/(m+1)) for j in 1:m,  k in 1:m]
	E, E2, -(√2/2)*V1, V
end

# ╔═╡ 77b036e6-7ddb-4a8e-8da5-3e0d1c9305ca
md"##  2.5 Block Multiplicationand Triangular Matrices"

# ╔═╡ 1a8d7932-4363-4749-ab31-a8ecf5288169
cm"""
4. If ``\boldsymbol{B}=\left[\boldsymbol{B}_1, \boldsymbol{B}_2\right]``, where ``\boldsymbol{B}_1 \in \mathbb{C}^{p \times r}`` and ``\boldsymbol{B}_2 \in \mathbb{C}^{p \times(n-r)}`` then
```math
\boldsymbol{A}\left[\boldsymbol{B}_1, \boldsymbol{B}_2\right]=\left[\boldsymbol{A} \boldsymbol{B}_1, \boldsymbol{A} \boldsymbol{B}_2\right]
```

5. If ``\boldsymbol{A}=\left[\begin{array}{l}\boldsymbol{A}_1 \\ \boldsymbol{A}_2\end{array}\right]``, where ``\boldsymbol{A}_1 \in \mathbb{C}^{k \times p}`` and ``\boldsymbol{A}_2 \in \mathbb{C}^{(m-k) \times p}`` then
```math
\left[\begin{array}{l}
\boldsymbol{A}_1 \\
\boldsymbol{A}_2
\end{array}\right] \boldsymbol{B}=\left[\begin{array}{l}
\boldsymbol{A}_1 \boldsymbol{B} \\
\boldsymbol{A}_2 \boldsymbol{B}
\end{array}\right]
```
"""

# ╔═╡ 4486fd72-e3f4-4e78-bc57-db8474dca974
cm"""
6. If ``\boldsymbol{A}=\left[\boldsymbol{A}_1, \boldsymbol{A}_2\right]`` and ``\boldsymbol{B}=\left[\begin{array}{l}\boldsymbol{B}_1 \\ \boldsymbol{B}_2\end{array}\right]``, where ``\boldsymbol{A}_1 \in \mathbb{C}^{m \times s}, \boldsymbol{A}_2 \in \mathbb{C}^{m \times(p-s)}, \boldsymbol{B}_1 \in`` ``\mathbb{C}^{s \times n}`` and ``\boldsymbol{B}_2 \in \mathbb{C}^{(p-s) \times n}`` then
```math
\left[\boldsymbol{A}_1, \boldsymbol{A}_2\right]\left[\begin{array}{l}
\boldsymbol{B}_1 \\
\boldsymbol{B}_2
\end{array}\right]=\left[\boldsymbol{A}_1 \boldsymbol{B}_1+\boldsymbol{A}_2 \boldsymbol{B}_2\right]
```
"""

# ╔═╡ 50ac8333-8f68-4f85-be14-dc99ba7fc398
cm"""
7. If ``\boldsymbol{A}=\left[\begin{array}{ll}\boldsymbol{A}_{11} & \boldsymbol{A}_{12} \\ \boldsymbol{A}_{21} & \boldsymbol{A}_{22}\end{array}\right]`` and ``\boldsymbol{B}=\left[\begin{array}{ll}\boldsymbol{B}_{11} & \boldsymbol{B}_{12} \\ \boldsymbol{B}_{21} & \boldsymbol{B}_{22}\end{array}\right]`` then
```math
\left[\begin{array}{ll}
\boldsymbol{A}_{11} & \boldsymbol{A}_{12} \\
\boldsymbol{A}_{21} & \boldsymbol{A}_{22}
\end{array}\right]\left[\begin{array}{ll}
\boldsymbol{B}_{11} & \boldsymbol{B}_{12} \\
\boldsymbol{B}_{21} & \boldsymbol{B}_{22}
\end{array}\right]=\left[\begin{array}{ll}
\boldsymbol{A}_{11} \boldsymbol{B}_{11}+\boldsymbol{A}_{12} \boldsymbol{B}_{21} & \boldsymbol{A}_{11} \boldsymbol{B}_{12}+\boldsymbol{A}_{12} \boldsymbol{B}_{22} \\
\boldsymbol{A}_{21} \boldsymbol{B}_{11}+\boldsymbol{A}_{22} \boldsymbol{B}_{21} & \boldsymbol{A}_{21} \boldsymbol{B}_{12}+\boldsymbol{A}_{22} \boldsymbol{B}_{22}
\end{array}\right],
```
"""

# ╔═╡ 519206e2-863d-4099-9f20-11cdad64a5e0
md"## 2.5.2 Triangular Matrices"

# ╔═╡ 0ec8c4e5-ba42-48c9-95ef-47457451e372
md"# Chapter 3: Gaussian Elimination and LU Factorizations"

# ╔═╡ adc9275a-493c-40c8-9b9a-05d38e43036e
cm"""
Thefollowing are called __elementary row operations__ applied on a  matrix:
1. Type 1: Exchanging two rows;    ``R_i \leftrightarrow`` ``R_j``
2. Type 2: Multiplying a row by a nonzero scalar; ``\lambda R_i \rightarrow R_i``
3. Type 3: Adding to a row a nonzero scalar multiple of another row. ``R_i+\lambda R_j \rightarrow`` ``R_i``.
"""

# ╔═╡ da8fd297-56f4-4175-8694-093de356dff9
md"## 3.2 Gauss and LU"

# ╔═╡ 71a38e11-f634-4e57-a277-e9eecd324381
cm"""
__Gaussian elimination without row interchanges__
(to slve a linear system ``\boldsymbol{A} \boldsymbol{x}=\boldsymbol{b}``) generates a sequence of equivalent systems 
- ``\boldsymbol{A}^{(k)} \boldsymbol{x}=\boldsymbol{b}^{(k)}`` for ``k=`` ``1, \ldots, n``, 

where ``\boldsymbol{A}^{(1)}=\boldsymbol{A}, \boldsymbol{b}^{(1)}=\boldsymbol{b}``, and ``\boldsymbol{A}^{(k)}`` has zeros under the diagonal in its first ``k-1`` columns. Thus ``\boldsymbol{A}^{(n)}`` is upper triangular and the system ``\boldsymbol{A}^{(n)} \boldsymbol{x}=\boldsymbol{b}^{(n)}`` is easy to solve. 

The matrix ``\boldsymbol{A}^{(k)}`` takes the form
```math
\boldsymbol{A}^{(k)}=\left[\begin{array}{ccc|ccccc}a_{1,1}^{(1)} & \cdots & a_{1, k-1}^{(1)} & a_{1, k}^{(1)} & \cdots & a_{1, j}^{(1)} & \cdots & a_{1, n}^{(1)} \\ & \ddots & \vdots & \vdots & & \vdots & & \vdots \\ & & a_{k-1, k-1}^{(k-1)} & a_{k-1, k}^{(k-1)} & \cdots & a_{k-1, j}^{(k-1)} & \cdots & a_{k-1, n}^{(k-1)} \\ \hline & & & a_{k, k}^{(k)} & \cdots & a_{k, j}^{(k)} & \cdots & a_{k, n}^{(k)} \\ & & & \vdots & & \vdots & & \vdots \\ & & & a_{i, k}^{(k)} & \cdots & a_{i, j}^{(k)} & \cdots & a_{i, n}^{(k)} \\ & & & \vdots & & \vdots & & \vdots \\ & & & a_{n, k}^{(k)} & \cdots & a_{n, j}^{(k)} & \cdots & a_{n, n}^{(k)}\end{array}\right].
```
"""

# ╔═╡ 8daf8b4c-39ed-47d2-932d-28b4c31ad6f7
let
	A = [1 1 0 3 
		 2 1 -1 1
		3 -1 -1 2
		-1 2 3 -1
	]
	b = [4;1;-3;4]
	L =[
		1 0 0 0
		2 1 0 0
		3 4 1 0
		-1 -3 0 1
	]
	U = [1 1 0 3;0 -1 -1 -5;0 0 3 13;0 0 0 -13]
	A , L*U
	# y =L\b
	# x = U\y
end

# ╔═╡ 420504bd-7c7b-4a09-a0c1-203b36ed340e
function lps(A,k)
	_,n = size(A)
	B = copy(A)
	for i in (k+1):n
		l=B[i,k]/B[k,k]
		for j in k:n
			B[i,j]=B[i,j]-l*B[k,j]
		end
	end
	B,B[1:k,1:k]
end

# ╔═╡ e9ab9e3f-47ae-48e7-ae3d-4afe1618919a
let
	A = [1 1 0 3 
		 2 1 -1 1
		3 -1 -1 2
		-1 2 3 -1
	]
	# b = [4;1;-3;4]
	# _,n =size(A) 
	# for k=1:1
	# 	for i in (k+1):n
	# 		l=A[i,k]/A[k,k]
	# 		for j in k:n
	# 			A[i,j]=A[i,j]-l*A[k,j]
	# 		end
	# 	end
	# end
	A1=A
	A2,=lps(A,1)
	A3, = lps(A2,2)
	A3
end

# ╔═╡ 6ed68d7b-7c50-4150-86f0-905fe6fd6a25
let
	A = [1 1 0 3 
		 2 1 -1 1
		3 -1 -1 2
		-1 2 3 -1
	]
	A1,Ak1=lps(A,1)
	A2,Ak2=lps(A1,2)
	A3,Ak3=lps(A2,3)
end

# ╔═╡ dff54fcc-2793-41d9-826a-8321e536b7d7
cm"""
- The complexity of an algorithm and say that the complexity of LU factorization of a full matrix is ``O\left(n^3\right)`` or more precisely ``\frac{2}{3} n^3``.
- that LU factorization is an ``O(n^3)``.
- solving a triangular system requires ``O(n^2)`` arithmetic operations.

"""

# ╔═╡ 953f4449-9dff-4275-b32b-075d727d2d2e
md"## 3.4 The PLU Factorization"

# ╔═╡ 2e68ed26-ae75-42f7-bc5d-fb23348e6711
md"### Pivoting"

# ╔═╡ d312c0d7-f80b-4da7-9ae3-9cdb2f783635
cm"""
- ``\boldsymbol{A} \boldsymbol{P}=\boldsymbol{A}(:, \boldsymbol{p})`` (Permutation of columns)
- ``\boldsymbol{P}^T \boldsymbol{A}=\boldsymbol{A}(\boldsymbol{p},:)`` (Permutation of rows)
- ``\boldsymbol{P}^T\boldsymbol{P} = \boldsymbol{I}=\boldsymbol{P}\boldsymbol{P}^T``

"""

# ╔═╡ 9521e54d-5a82-44f8-a4b6-94b229a58ec5
let
	Random.seed!(123)
	A = rand(-3:8,3,3)
	perms = collect(permutations(1:3))
	p = rand(perms,1)[1]
	A, A[:,p], A[p,:]
	P=I(3)[:,p]
	# P=P[p,:]
	# P*P'
end

# ╔═╡ b9a41148-048e-463c-b0c1-85de79ec6ee1
cm"""
- ``\boldsymbol{I}_{i k}^2=\boldsymbol{I}`` (interchange matrix is symmetric and equal to its own inverse.)
- We can keep track of the row interchanges using pivot vectors ``\boldsymbol{p}_k``. We define 
```math
\boldsymbol{p}:=\boldsymbol{p}_n, \text{ where } \boldsymbol{p}_1:=[1,2, \ldots, n]^T, \text{ and } \boldsymbol{p}_{k+1}:=\boldsymbol{I}_{r_k, k} \boldsymbol{p}_k \text{ for } k=1, \ldots, n-1
```
"""

# ╔═╡ aa5e1849-e9e4-4b2f-90b7-9716ad441a85
cm"""
```math
\boldsymbol{P}^T=\boldsymbol{P}_{n-1} \cdots \boldsymbol{P}_1=\boldsymbol{I}(\boldsymbol{p},:), \quad \boldsymbol{P}=\boldsymbol{P}_1 \boldsymbol{P}_2 \cdots \boldsymbol{P}_{n-1}=\boldsymbol{I}(:, \boldsymbol{p})
```

Instead of interchanging the rows of ``\boldsymbol{A}`` during elimination we can keep track of the ordering of the rows using the pivot vectors ``\boldsymbol{p}_k``. Gaussian elimination with row pivoting starting with ``a_{i j}^{(1)}=a_{i j}`` can be described as follows:
```math
\begin{aligned}
& \boldsymbol{p}=[1, \ldots, n]^T \\
& \text { for } k=1: n-1 \\
& \text { choose } r_k \geq k \text { so that } a_{p_{r_k}, k}^{(k)} \neq 0 \\
& \quad \boldsymbol{p}=I_{r_k, k} \boldsymbol{p} \\
& \text { for } i=k+1: n \\
& \qquad a_{p_i, k}^{(k)}=a_{p_i, k}^{(k)} / a_{p_k, k}^{(k)} \\
& \text { for } j=k: n \\
& \qquad a_{p_i, j}^{(k+1)}=a_{p_i, j}^{(k)}-a_{p_i, k}^{(k)} a_{p_k, j}^{(k)}
\end{aligned}
```
"""

# ╔═╡ 8816f81b-3f55-49e0-ae40-ff6672c951a9
cm"""
- The PLU factorization can also be written ``\boldsymbol{P}^T \boldsymbol{A}=\boldsymbol{L} \boldsymbol{U}``.
"""

# ╔═╡ 2d53670e-772e-41b3-9076-3193b40a5c09
md"###  Pivot Strategie"

# ╔═╡ 12e05a32-7a96-420c-a90b-1b76748cb09f
cm"""
In __partial pivoting__ we select the __largest element__
```math
\left|a_{r_k, k}^{(k)}\right|:=\max \left\{\left|a_{i, k}^{(k)}\right|: k \leq i \leq n\right\}
```
with ``r_k`` the smallest such index in case of a tie. The following example illustrating that small pivots should be avoided.
"""

# ╔═╡ b96f87dd-14b9-49ce-970a-b5093fc13e94
# let
# 	A = [1 1 0 3 
# 		 2 1 -1 1
# 		3 -1 -1 2
# 		-1 2 3 -1
# 	]
# 	F = lu(A)
# 	L1,U1,pt = F.L, F.U, F.p 
# 	P1 = I(4)[:,pt]
# 	A≈ P1*L1*U1
# 	B = copy(Rational.(A))
# 	B[1,:],B[3,:]=B[3,:],B[1,:]
# 	B[2,:]=-(2//3)*B[1,:]+B[2,:]
# 	B[3,:]=-(1//3)*B[1,:]+B[3,:]
# 	B[4,:]=(1//3)*B[1,:]+B[4,:]
# 	B[3,:]=-(4//5)*B[2,:]+B[3,:]
# 	B[4,:]=-B[2,:]+B[4,:]
# 	B[3,:],B[4,:]=B[4,:],B[3,:]
# 	B[4,:]=-1//5*B[3,:]+B[4,:]
# 	B
# 	U = B
# 	L =[1 0 0 0;
# 		2//3 1 0 0;
# 		-1//3 1 1 0;
# 		1//3 4//5 1//5 1]
# 	p2=[4,2,1,3]
# 	P2 = I(4)[p2,:]
# 	A,P2*L*B
# 	# B,U1
# 	P2
# 	# p2,pt
# 	# U1,Float64.(U)
	
# end

# ╔═╡ a78b0614-3a52-4fb8-acc7-147707e1d456
let
	A=[1 -1 1 2
	-2 1 1 1
	2 -1 2 3
	-4 1 0 2
	]
	B = Rational.(copy(A))
	B[4,:],B[1,:]=B[1,:],B[4,:]
	B[2,:] = -(1//2)*B[1,:]+B[2,:]
	B[3,:] = (1//2)*B[1,:]+B[3,:]
	B[4,:] = (1//4)*B[1,:]+B[4,:]
	B[2,:],B[4,:]=B[4,:],B[2,:]
	B[3,:] = -(2//3)*B[2,:]+B[3,:]
	B[4,:] = (2//3)*B[2,:]+B[4,:]
	B[3,:],B[4,:]=B[4,:],B[3,:]
	B[4,:] = -(4//5)*B[3,:]+B[4,:]
	pt =[4,1,2,3]
	P = I(4)[:,pt]
	# L =[
	# 	1 0 0 0
	# 	-1//4 1 0 0
	# 	 1//2 -2//3 1 0
	# 	-1//2 2//3 4//5 1
	# ]
	# P*L*B , A
	L,U,p = lu(Rational.(A),NoPivot())
	# P=I(4)[p,:]
	# L=L*diagm(0=>diag(U))
	
	# U[1,:]=U[1,:]//U[1,1]
	# U[2,:]=U[2,:]//U[2,2]
	# U[3,:]=U[3,:]//U[3,3]
	# U[4,:]=U[4,:]//U[4,4]
	
	# L*U,A
	L, U
	# p
	# P=I(4)[p,:]
	# p
	# B
end

# ╔═╡ fd98de8b-bb26-4c7c-8b65-d16f3fe57a09
cm"""
Use Julia `lu` for PLU factorization:
```julia

L,U,p = lu(A)
P = I[p,:]
```
"""

# ╔═╡ 7f2f4c83-63c9-401f-b6c5-25722674a68c
md"## 3.5 The LU and LDU Factorizations"

# ╔═╡ c4ec32de-d6e1-4ef4-b870-cd5c5e7f77fa
cm"""
Gaussian elimination without row interchanges is one way of computing an LU factorization of a matrix. 
There are other ways that can be advantageous for certain problems:
- __L1U__: ``\quad l_{i i}=1`` all ``i``,
- __LU1__: ``u_{i i}=1`` all ``i``,
- __LDU__: ``\quad \boldsymbol{A}=\boldsymbol{L} \boldsymbol{D} \boldsymbol{U}, l_{i i}=u_{i i}=1`` all ``i, \boldsymbol{D}=\operatorname{diag}\left(d_{11}, \ldots, d_{n n}\right)``.
"""

# ╔═╡ d5dd27d9-16dc-44b0-9128-63e16cc9db10
let
	Random.seed!(286)
	A = rand(-2:0.2:4,10,10)
	j =rand(1:10,10)
	foreach(j) do i
	A[i,i]=0
	end
	A
	minors =[det(A[1:i,1:i]) for i in 1:9]
end

# ╔═╡ 93c8a52c-e3d9-42a5-b9be-11e6a090fa55
md"## 3.6 Block LU Factorization"

# ╔═╡ c9344f33-0372-48de-ae18-819e05d5852f
cm"""
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is a block matrix of the form
```math
\boldsymbol{A}:=\left[\begin{array}{ccc}
\boldsymbol{A}_{11} & \cdots & \boldsymbol{A}_{1 m} \\
\vdots & & \vdots \\
\boldsymbol{A}_{m 1} & \cdots & \boldsymbol{A}_{m m}
\end{array}\right]
```
where each diagonal block ``\boldsymbol{A}_{i i}`` is square. We call the factorization
```math
\boldsymbol{A}=\boldsymbol{L} \boldsymbol{U}=\left[\begin{array}{ccccc}
\boldsymbol{I} & & & \\
\boldsymbol{L}_{21} & \boldsymbol{I} & & \\
\vdots & & \ddots & \\
\boldsymbol{L}_{m 1} & \cdots & \boldsymbol{L}_{m, m-1} & \boldsymbol{I}
\end{array}\right]\left[\begin{array}{cccc}
\boldsymbol{U}_{11} & & \cdots & \boldsymbol{U}_{1 m} \\
& \boldsymbol{U}_{22} & \cdots & \boldsymbol{U}_{2 m} \\
& & \ddots & \vdots \\
& & \boldsymbol{U}_{m m}
\end{array}\right]
```
"""

# ╔═╡ 91739d93-dfda-4599-9c4b-68c6ee7f8df1
md"# Chapter 4: LDL* Factorization and Positive Definite Matrices"

# ╔═╡ e440727e-955e-4db6-a9d5-c1125917f56e
md"##  4.1 The LDL* sFactorization"

# ╔═╡ 43ffc137-855c-4e16-8ebd-ce2f73ad7161
md"## 4.2 Positive Definite and Semidefinite Matrices"

# ╔═╡ aee4e1e3-7ce8-48a9-b258-b0ab46d52925
let
	Random.seed!(123)
	x = rand(Complex{Int8},2,1)
	y = rand(Complex{Int8},2,1)
	conj(x⋅y),y⋅x
end

# ╔═╡ 380db7d1-9bb6-44e3-9063-6d98bceaf15f
cm"""
- If ``f:\Omega\subset \mathbb{R}^n\to \mathbb{R}``, then
```math
\nabla f(\boldsymbol{x})=\left[\begin{array}{c}\frac{\partial f(\boldsymbol{x})}{\partial x_1} \\ \vdots \\ \frac{\partial f(\boldsymbol{x})}{\partial x_n}\end{array}\right] \in \mathbb{R}^n, \quad H f(\boldsymbol{x})=\left[\begin{array}{ccc}\frac{\partial^2 f(\boldsymbol{x})}{\partial x_1 \partial x_1} & \cdots & \frac{\partial^2 f(\boldsymbol{x})}{\partial x_n \partial x_1} \\ \vdots & \vdots \\ \frac{\partial^2 f(\boldsymbol{x})}{\partial x_1 \partial x_n} & \ldots & \frac{\partial^2 f(\boldsymbol{x})}{\partial x_n \partial x_n}\end{array}\right] \in \mathbb{R}^{n \times n}
```
"""

# ╔═╡ faf07e34-58fd-4be3-a81d-9fddd668e582
md"##  4.2.1 The Cholesky Factorization"

# ╔═╡ e9ea8099-acfd-4df4-a19a-f510888e1f98
let
	
	Ao = Rational.([2 4 -3;4 14 -9;-3 -9 12])
	A = deepcopy(Ao)
	# A[2,:] = -2A[1,:]+A[2,:]
	# A[3,:] = (3//2)A[1,:]+A[3,:]
	# A[3,:] = (1//2)A[2,:]+A[3,:]
	A
	# L1,U1 = lu(Ao,NoPivot())
	# L1,U1
	# D =diagm(0=>[2;6;6])
	# U=L1'
	# A = L1*D*U
	# L,U = cholesky(Ao)
	# L,U
	# L, U
	# U12=deepcopy(U1)
	# D = diagm(0=>diag(U1))
	# U12=stack(map(x->U1[:,x[1]] ./ diag(U1),enumerate(eachrow(U1))))
	# L1,U12
	
	# U12',L1
end

# ╔═╡ 0275d178-b918-445a-be9e-eac9c7aaa471
md"## 4.2.2 Positive Definite and Positive Semidefinite Criteria"

# ╔═╡ 6d988a2c-f6ae-418d-9cd4-11207ec7a2ce
cm"""
Consider 
```math
A= \begin{bmatrix}
20 & 6 - 3i & 1 + 1i & 4 + 4i \\ 
6 + 3i & 24 & 6 - 1i & 6 + 1i \\ 
1 - 1i & 6 + 1i & 20 & -4 - 2i \\ 
4 - 4i & 6 - 1i & -4 + 2i & 16  
\end{bmatrix}
```
"""

# ╔═╡ 02c7ec2d-e57f-49d5-a8c3-a6f0637e026f
let
	Random.seed!(122)
	# A = rand(-2:4,4,4)
	A = [20 + 0im 6 - 3im 1 + 1im 4 + 4im; 6 + 3im 24 + 0im 6 - 1im 6 + 1im; 1 - 1im 6 + 1im 20 + 0im -4 - 2im; 4 - 4im 6 - 1im -4 + 2im 16 + 0im]
	
	# isposdef(A)
	# sqrt(24*20),4*sqrt(30),3sqrt(5)
end

# ╔═╡ 41526c5c-0be0-4bef-af4c-557c28b918cd
md"# Chapter 5: Orthonormal and Unitary Transformations"
 # 5.1− 5.4, 5.6

# ╔═╡ 694c5383-28cf-4e69-90bc-366a40c1e0a9
md"##  5.1 Inner Products, Orthogonality and Unitary Matrices"

# ╔═╡ b012a15b-be75-4cb4-96ff-d679ed4376d1
md"### Orthogonality"

# ╔═╡ 8509b4b9-b4c8-47e5-84ae-c21b47eafa16
let
	s1 =[1;2;-2]
	s2=[1;3;2]
	s3 =[2;5;1]
	# det([s1 s2 s3])
	s3 	⋅s2
# 	v1=s1
# 	v2 = s2 -v1*(s2⋅v1)/(norm(v1))^2
# 	v3 = s3 - v1*(s3⋅v1)/(norm(v1))^2 - v2*(s3⋅v2)/(norm(v2))^2
# 	v1=normalize(v1)
# 	v2=normalize(v2)
# 	v3=normalize(v3)
# 	norm.(eachcol([v1 v2 v3]))
end

# ╔═╡ 2bebb2e6-55b6-4a33-9409-2b901b6dfc5e
md"### Sum of Subspaces and Orthogonal Projections"

# ╔═╡ 5e8ea837-cd68-4dd1-9d4b-925223fc63f9
cm"""
Suppose ``\mathcal{S}`` and ``\mathcal{T}`` are subspaces of a real or complex vector space ``\mathcal{V}`` endowed with an inner product ``\langle\boldsymbol{x}, \boldsymbol{y}\rangle``. We define
- Sum: ``\mathcal{S}+\mathcal{T}:=\{s+t: s \in \mathcal{S}`` and ``t \in \mathcal{T}\}``,
- ``\operatorname{direct} \operatorname{sum} \mathcal{S} \oplus \mathcal{T}`` : a sum where ``\mathcal{S} \cap \mathcal{T}=\{0\}``,
- orthogonal sum ``\mathcal{S} \stackrel{\perp}{\oplus} \mathcal{T}`` : a sum where ``\langle\boldsymbol{s}, \boldsymbol{t}\rangle=0`` for all ``s \in \mathcal{S}`` and ``\boldsymbol{t} \in \mathcal{T}``.
"""

# ╔═╡ 17312234-ed45-4e26-868e-3f25c14f73bd
cm"### Unitary and Orthogonal Matrices"

# ╔═╡ 3d678fbd-65f7-4516-a684-6724c66970de


# ╔═╡ a916bede-5304-4120-a929-979d8fbff63c
md"## 5.2 The Householder Transformation"

# ╔═╡ b5ba03bb-8c78-4a14-91d8-f4b110db8885


# ╔═╡ 9d950ca7-0162-41a6-b59f-ba2d64d48f4e
cm"""
- If ``\boldsymbol{x}=0`` then ``\boldsymbol{H} \boldsymbol{x}=\mathbf{0}`` and ``a=0``. Any ``\boldsymbol{u}`` with ``\|\boldsymbol{u}\|_2=\sqrt{2}`` will work, and we choose ``\boldsymbol{u}:=\sqrt{2} \boldsymbol{e}_1`` in this case. 
- For ``\boldsymbol{x} \neq \mathbf{0}`` we define
```math
\boldsymbol{u}:=\frac{\boldsymbol{z}+\boldsymbol{e}_1}{\sqrt{1+z_1}}, \text { where } \boldsymbol{z}:=\bar{\rho} \boldsymbol{x} /\|\boldsymbol{x}\|_2 .
```
"""

# ╔═╡ 7e6ce1a2-8de1-4170-9585-62b25311e646
begin
	function housegen(x::Vector{T}) where T
	    # Initialize variables
	    a = norm(x)
	    
	    # If the norm is zero, handle this edge case
	    if a == 0
	        u = copy(x)
	        u[1] = sqrt(2)
	        return u, a
	    end
	
	    # Determine r based on the sign of x[1]
	    if x[1] == 0
	        ρ = one(T)  # ρ = 1
	    else
	        ρ = x[1] / abs(x[1])
	    end
	
	    # Compute u and a
	    u = conj(ρ) * x / a
	    u[1] = u[1] + 1
	    u = u/sqrt(abs(u[1]))
	    a = -ρ * a
	    
	    return u, a
	end
	H(u)= x->x-(u'*x)*u
end

# ╔═╡ d94ca4da-f011-4c8c-b92e-985d29c4f3e5
let
	x =[-1,2]
	u,a = housegen(x)
	
	H(u)(x)
end

# ╔═╡ 0865697a-e057-4be3-9aa9-a1bf8831829d
cm"""
 - Householder transformations can also be used to zero out only the lower part of
 a vector.
"""

# ╔═╡ 7c3661da-0ee8-4877-8d50-07c6704dbbb0
let
	y=[1,2]
	z = [-2,3]
	x = [y;z]
	us,a = housegen(z)
	u=[zeros(2);us]
	u'*u
	H(u)(x),	a*[1,0]
end

# ╔═╡ 8acef7cd-a100-46be-b3ce-6b00ab56f479
md"## 5.3 Householder Triangulation"

# ╔═╡ 864a50be-eabc-44be-891e-64d7caa2d8ef
md"### The Algorithm"

# ╔═╡ 7fd8fae9-c083-46c7-a0b9-1c8317641669
cm"""
Let ``A\in \mathbb{C}^{m\times n}``. We need to transform it into upper trapezoidal form using Householder transformations.

- If ``m>n`` (TALL MATRIX), We find 
```math
\boldsymbol{A}_{n+1}:=\boldsymbol{H}_n \boldsymbol{H}_{n-1} \cdots \boldsymbol{H}_1 \boldsymbol{A}=\left[\begin{array}{c}
\boldsymbol{R}_1 \\
\mathbf{0}
\end{array}\right]=\boldsymbol{R}
```
and where ``\boldsymbol{R}_1`` is square and upper triangular. We define
```math
\boldsymbol{A}_1:=\boldsymbol{A}, \quad \boldsymbol{A}_{k+1}=\boldsymbol{H}_k \boldsymbol{A}_k, \quad k=1,2, \ldots, n
```
"""

# ╔═╡ b35a0a84-880b-4da3-a01d-f7d7fd1bb2a4
cm"""
Suppose ``\boldsymbol{A}_k`` has the following form
```math
\begin{array}{lclc}
\boldsymbol{A}_k&=&\left[\begin{array}{ccc|ccccc}a_{1,1}^{(1)} & \cdots & a_{1, k-1}^{(1)} & a_{1, k}^{(1)} & \cdots & a_{1, j}^{(1)} & \cdots & a_{1, n}^{(1)} \\ & \ddots & \vdots & \vdots & & \vdots & & \vdots \\ & & a_{k-1, k-1}^{(k-1)} & a_{k-1, k}^{(k-1)} & \cdots & a_{k-1, j}^{(k-1)} & \cdots & a_{k-1, n}^{(k-1)} \\ \hline & & & a_{k, k}^{(k)} & \cdots & a_{k, j}^{(k)} & \cdots & a_{k, n}^{(k)} \\ & & & \vdots & & \vdots & & \vdots \\ & & & a_{i, k}^{(k)} & \cdots & a_{i, j}^{(k)} & \cdots & a_{i, n}^{(k)} \\ & & & \vdots & & \vdots & & \vdots \\ & & & a_{m, k}^{(k)} & \cdots & a_{m, j}^{(k)} & \cdots & a_{m, n}^{(k)}\end{array}\right]\\
&=&\left[\begin{array}{cc}\boldsymbol{B}_k & \boldsymbol{C}_k \\ \mathbf{0} & \boldsymbol{D}_k\end{array}\right]
\end{array}
```

"""

# ╔═╡ 2d60f5a5-139f-4451-af2e-33ee04355709
cm"""
__Input__: Matrices ``A`` and ``B``

__Output__: Matrices ``R`` and ``C``
    
__Step 1__: Get the sizes of ``A`` and ``B``
	- m = number of rows in ``A``
    - n = number of columns in ``A``
    - r = number of columns in ``B``

__Step 2__: Concatenate ``A`` and ``B`` horizontally into a new matrix

__Step 3__: Perform Householder transformation
        
For each ``k`` from 1 to ``\min(n, m-1)``:
- a. Compute the Householder vector ``v`` for the submatrix ``A(k:m, k)``
- b. Update the submatrix ``A(k:m, k:n+r)`` using the Householder reflection

__Step 4__: Extract the upper triangular part of ``A`` and assign it to ``R``

__Step 5__: Extract the right part of ``A`` (from column ``n+1`` to ``n+r``) and assign it to ``C``
    
"""

# ╔═╡ dca4ccc0-bb15-41b3-9fb5-a3d74bf552c9
begin
	function housetriang(A)
		m, = size(A)
		housetriang(A,I(m))
	end
	function housetriang(A,B)
		m,n = size(A)
		r=size(B,2)
		A=[A B]
		for k in 1:min(n,m-1)
			v, a = housegen(A[k:m,k])
			A[k,k]=a
			C = A[k:m,k+1:n+r]
			A[k:m,k+1:n+r]=C-v*(v'*C)
		end
		R = triu(A[:,1:n])
		C = A[:,n+1:n+r]
		R,C
	end
end

# ╔═╡ f79e868a-5ed7-4439-ac3c-9229f64360d2
let
	A = rand(-1:6,6,4)
	R, C = housetriang(Float64.(A))
	R
end

# ╔═╡ 04b6d34a-7a30-4337-a862-a84f17680ac6
md"### Solving Linear Systems Using Unitary Transformations"

# ╔═╡ 3ba9e235-869c-4011-afeb-c09385260140
let
	A = [1.0  2  1  
	3  8  7     
	2  7  9  
	]
	b= [4;20;23.0]
	R, C = housetriang(A)
	C
	bb= C*b
	R
	x3 = bb[end]/R[end,end]
	x2 = (bb[2]-R[2,3]*x3)/R[2,2]
	x1 = (bb[1]-R[1,3]*x3-R[1,2]*x2)/R[1,1]
	A*[x1;x2;x3]-b
	[x1,x2,x3],A\b
end

# ╔═╡ 3e8a3007-5161-40f5-bf7d-47d8aeecb1aa
let
	Au = [3 1 -1  14 
	 1 -2 5   -7 
	 4  1 2  17.0]
	R,C = housetriang(Au)
end

# ╔═╡ 7ea5927a-81bb-404a-aa62-67bce4e00313
md"##  5.4 The QR Decomposition and QR Factorization"

# ╔═╡ 211f7f09-31f4-44d5-8e8f-6bb5a2920ae7
cm"""
- Gaussian elimination without row interchanges results in an LU factorization 

```math
\boldsymbol{A}=\boldsymbol{L} \boldsymbol{U}\text{ of }\boldsymbol{A} \in \mathbb{C}^{n \times n}
``` 
- Householder triangulation of ``\boldsymbol{A}``. Applying Algorithm above gives 

```math
\boldsymbol{R}=\boldsymbol{H}_{n-1} \cdots \boldsymbol{H}_1 \boldsymbol{A}
```
- That is 
```math
\boldsymbol{A}=\boldsymbol{Q} \boldsymbol{R}, \text{ where } \boldsymbol{Q}=\boldsymbol{H}_1 \cdots \boldsymbol{H}_{n-1}
\text{ is unitary and }\boldsymbol{R} \text{ is upper triangular.}
```

> This is known as a __QR factorization__ of ``\boldsymbol{A}``.
"""

# ╔═╡ ca414318-5a80-423e-b4a3-561e282e4111
let
	A=[1  3  1. 
	1  3  7 
	1  -1  -4 
	1  -1  2]
	Q,R = qr(A)
	collect(Q)
	# Q*R
	Q[:,1:3]*R[1:3,:]
	# R
end

# ╔═╡ 2bb49eb8-1c93-44f3-a05b-745594830356
# solve


# ╔═╡ cc72a9e5-1f40-496e-8f87-79cf6734da89
md"##  5.5 QR and Gram-Schmid"

# ╔═╡ fbf9565f-f1d6-4a50-b81d-4863eff22d3b
begin
	# function QRFactor(A)
	# 	m,n =size(A)
	# 	vjs = Vector{Vector{<:Real}}(undef,n)
	# 	foreach(enumerate(eachcol(A))) do (j,a)
	# 		 vj=if j==1
	# 			 a
	# 		 else
	# 			 a - sum(vjs[i]*dot(vjs[i],a)/dot(vjs[i],vjs[i]) for i in 1:j-1)
	# 		 end
	# 		vjs[j]=[vj...]
	# 	end
	# 	vjs_norms = map(norm,vjs)
	# 	Q1 = stack(vjs ./ vjs_norms)
	# 	dict = [i=>map(j->dot(A[:,j+i],Q1[:,j]),1:n-i) for i in 1:n-1]
	# 	R1 =diagm(0=>vjs_norms, dict...)
		
	# 	Q1, R1
	# end
	function QRFactor(A)
		m,n = size(A)
		Q = copy(A)
		R = zeros(eltype(A),n,n)
		for k in 1:n
 			R[1:k-1,k] = Q[:,1:k-1]'*A[:,k]
 			Q[:,k] = A[:,k]- Q[:,1:k-1]*R[1:k-1,k]
 			R[k,k] = norm(Q[:,k])
 			Q[:,k] = Q[:,k]/R[k,k]
 		end
		Q,R
	end

end

# ╔═╡ c6576e1e-f86f-46a4-bd2f-7b82550c231d
let
	A=[1  3  1. 
	1  3  7 
	1  -1  -4 
	1  -1  2]
	Q1,R1= QRFactor(A)
	Q,R = qr(A)
	collect(Q),Q1,R,R1
	# dot(Q[:,2],Q[:,3])
	Q1*R1

end

# ╔═╡ 22cac4d7-7c2a-43e9-8d94-ed6df2508fed
md"# Chapter 6: Eigenpairs and Similarity Transformations"

# ╔═╡ 5a573a95-6092-4def-b900-1b80c3aca31f
md"## 6.1 Defective and Nondefective Matrices"

# ╔═╡ e7a24ffb-1a85-433c-a9c2-4868f9ba4205
md"### Similarity Transformations"

# ╔═╡ 2aa4b07a-f768-4986-814d-eb5c6a279adc
cm"""
1. Similar matrices have the same eigenvalues, they even have the same charac
teristic polynomial.
2. ( ``\boldsymbol{\lambda , \boldsymbol { x }}`` ) is an eigenpair for ``\boldsymbol{S}^{-1} \boldsymbol{A} \boldsymbol{S}`` if and only if ``(\boldsymbol{\lambda}, \boldsymbol{S} \boldsymbol{x})`` is an eigenpair for ``\boldsymbol{A}``
3. If ``\boldsymbol{S}^{-1} \boldsymbol{A} \boldsymbol{S}=\boldsymbol{D}=\operatorname{diag}\left(\lambda_1, \ldots, \lambda_n\right)`` then the columns of ``\boldsymbol{S}`` are eigenvectors of ``\boldsymbol{A}`` and ``\boldsymbol{A}`` is nondefective . 
<div style="padding-left:40px;">

Conversely, if ``\boldsymbol{A}`` is nondefective then it can be diagonalized by a similarity transformation ``\boldsymbol{S}^{-1} \boldsymbol{A}\boldsymbol{S}``, where the columns of ``\boldsymbol{S}`` are eigenvectors of ``\boldsymbol{A}``.
</div>
"""

# ╔═╡ bf91a78d-429c-46c1-8a09-8f7551c32d6a
md"### Algebraic and Geometric Multiplicity of Eigenvalues"

# ╔═╡ 1660650f-917c-4354-acb7-c287f53548b8
let
	A = [6  3  4 
	0  6  2 
	0  0  7]
	λs, Λ = eigen(A)
	# Λ[:,3], normalize([10,2,1])

	
	 
end

# ╔═╡ 8221aa32-82ae-4ec6-b40b-0ec7c54f92de
cm"""
- ``A`` is not diagnalizable
- ``A`` is defective
"""

# ╔═╡ 3ea9b542-6182-48d3-9fc6-29b7d1c28920
md"## 6.2 The Jordan Factorization"

# ╔═╡ 43dc89bb-3c51-4aea-98ee-4df58ccd2997
let
	A = [6  3  4 
	0  6  2 
	0  0  7]

	N1 =A-6I
	# N1^2
	v1 =[3;0;0]
	v2 =[0;1;0]
	v3=[10;2;1]
	S = [v1 v2 v3]
	J =[6 1 0
		0 6 0
		0 0 7
	]	
	S*J*inv(S), A
	
end

# ╔═╡ df706f98-7cf7-4b84-9868-0e67e9e4c34a
let
	A=[2  1  -1  8  -3 
	0  2  0  7  5 
	0  0  2  7  5 
	0  0  0  2  0 
	0  0  0  0  2.0
	]
	N = A-2I
	v1=[1;0;0;0;0]
	v2=[0;1;1;0;0]
	v3=[0;61;0;-5;7]
	N^2
	w3 = [0;1;0;0;0]
	N*w3
	w4 = [0;0;0;1;0]
	v2= N*w4
	S =[v1 w3 v2 w4 v3]
	J =[
		2 1 0 0 0
		0 2 0 0 0
		0 0 2 1 0
		0 0 0 2 0
		0 0 0 0 2
	]
	round.(S*J*inv(S)),A
end

# ╔═╡ 6e123b71-f3bd-40b5-ad65-82515b0a44df
let
	# # @variables λ::Real
	# A=[2  0  0  0 
	# 0  2  0  0 
	# 0  0  2  1 
	# 1  0  0  2]
	# # det(A-λ*I)
	# N=A-2I
	# N^3
	# # level 1 N
	# v1=[0;1;0;0]
	# v2=[0;0;1;0]
	# # level 2 N^2
	# v3=[0;0;0;1]

	# # level 3 N^3
	# v4=[1;0;0;0]
	# N*v4
	# N^2*v4
	# S= [v2 v3 v4 v1]
	# J =[2 1 0 0
	# 	0 2 1 0
	# 	0 0 2 0
	# 	0 0 0 2]
	# S*J*inv(S)
	
end

# ╔═╡ 63c34957-6379-4723-b537-f9753f6fe814
md"## 6.3 The Schur Factorization and Normal Matrices"

# ╔═╡ 6be7cbc4-023b-4d2e-bda5-36a996fad82c
md"# Chapter 7: The Singular Value Decomposition (svd)"

# ╔═╡ 2f9a1b80-5247-4af9-bf89-15a50eb7e5c6
md"##  Normal Matrices"

# ╔═╡ 4bf58173-2848-4256-bc2d-01f0e3136094
md"## 7.1 The SVD Always Exists
### The Matrices $A^* A, A A^*$"

# ╔═╡ 5d9e665a-7dd6-42cb-99f8-47ab69256feb
let
	A = [1 1 0.0;0 0 1]
	# B = A'*A
	# eigen(B)
	# w1 =[1;1;0]
	# v1=w1/norm(w1)
	# w2 =[0;0;1]
	# v2=w2/norm(w2)
	# w3 =[-1;1;0]
	# v3=w3/norm(w3)
	# V =[v1 v2 v3]
	# u1,u2 = A*v1/sqrt(2),A*v2
	# U = [u1 u2]
	# Σ = [sqrt(2) 0 0;0 1 0]
	# U*Σ*V', A
	# V*Σ'U'
	
	
end

# ╔═╡ f4cf76b0-e5f6-4920-93a7-0be1f1041987
let
	A = [1  0  1 
		 0  1  0 
		 0  1  1 
		 0  1  0 
		 1  1  0]
	B = A'*A
	# @variables λ::Real
	# expand(det(B-λ*I))
	# λ = 5,2,1
	# L,U=lu(B-λ[1]I,NoPivot())
	# U
	
	# w1=[1;2;1]
	# v1=w1/norm(w1)

	# B-λ[2]I
	# w2=[1;-1;1]
	# v2=w2/norm(w2)

	# B-λ[3]I
	# w3=[-1;0;1]
	# v3=w3/norm(w3)
	
	# u1=(1/√λ[1])A*v1
	# norm(u1)
	# u2 = (1/√λ[2])*A*v2
	# # norm(u2)
	# u3 = sqrt(2)*A*v3
	# # u3=u3/norm(u3)
	# # C = A*A'
	
	# u4 =[0;-1;0;1;0]
	# u4 =u4/norm(u4)
	# u5 =[-1;-2;1;0;1]
	# u5=u5-dot(u4,u5)*u4
	# u5=u5/norm(u5)
	# U =[u1 u2 u3 u4 u5]
	# V = [v1 v2 v3]
	# S=[sqrt(5) 0 0
	#    0 sqrt(2) 0
	#    0 0 1
	# 	0 0 0
	# 	0 0 0
	# ]
	# round.(U*S*V')
	# U2,S2,Vt2 = svd(A)
	# diag(S),S2
	# dot(u4,u5)
	
end

# ╔═╡ 57322645-33db-4c87-a5c8-51d69e8e74e9
cm"""

__The Singular Value Factorization__

Suppose ``\boldsymbol{A}=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*`` is a singular value decomposition of ``\boldsymbol{A}`` of rank ``r``. Consider the block partitions
```math
\begin{aligned}
& \boldsymbol{U}=\left[\boldsymbol{U}_1, \boldsymbol{U}_2\right] \in \mathbb{C}^{m \times m}, \quad \boldsymbol{U}_1:=\left[\boldsymbol{u}_1, \ldots, \boldsymbol{u}_r\right], \quad \boldsymbol{U}_2:=\left[\boldsymbol{u}_{r+1}, \ldots, \boldsymbol{u}_m\right], \\
& \boldsymbol{V}=\left[\boldsymbol{V}_1, \boldsymbol{V}_2\right] \in \mathbb{C}^{n \times n}, \quad \boldsymbol{V}_1:=\left[\boldsymbol{v}_1, \ldots, \boldsymbol{v}_r\right], \quad \boldsymbol{V}_2:=\left[\boldsymbol{v}_{r+1}, \ldots, \boldsymbol{v}_n\right], \\
& \boldsymbol{\Sigma}=\left[\begin{array}{cc}
\boldsymbol{\Sigma}_1 & \mathbf{0}_{r, n-r} \\
\mathbf{0}_{m-r, r} & \mathbf{0}_{m-r, n-r}
\end{array}\right] \in \mathbb{R}^{m \times n}, \text { where } \boldsymbol{\Sigma}_1:=\operatorname{diag}\left(\sigma_1, \ldots, \sigma_r\right)
\end{aligned}
```
So
```math
\boldsymbol{A}=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*=\boldsymbol{U}_1 \boldsymbol{\Sigma}_1 \boldsymbol{V}_1^*.
```
"""

# ╔═╡ 4843fc90-c81f-4810-8da8-8b1ce0694a55
md"## SVD and the Four Fundamental Subspaces"

# ╔═╡ 2a8a27d1-55ff-4471-80c3-28a7c605b97e
md"## 7.3 A Geometric Interpretation"

# ╔═╡ 394b491d-36b0-48a8-92d0-b2abb1ffdc01
let
	A= [11  48 ; 48  39]/25
	U,S,Vt = svd(A)
	
	a = S[1]     # Semi-major axis
	b = S[2]    # Semi-minor axis
	cx, cy = 0, 0  # Center of the ellipse
	θ = π/4   # Rotation angle (in radians)
	t = range(0, 2π, length=500)
	x = a * cos.(t)
	y = b * sin.(t)
	
	# Rotate the ellipse by θ
	x_rot = U * [x';y']
	# y_rot = sin(θ) * x + cos(θ) * y
	
	# Shift to the desired center
	# x_shifted = x_rot .+ cx
	# y_shifted = y_rot .+ cy
	
	# Plot the filled ellipse
	p1 = plot(x, y, fill = true, lw=2, label="", color=:lightblue, alpha=0.5)
	plot(p1,frame_style=:origin)
	p2 = plot(p1,x_rot[1,:], x_rot[2,:], fill = true, lw=2, label="", color=:lightgreen, alpha=0.5)
	plot(p2,frame_style=:origin)
	# plot(p1,p2)
	
	
end

# ╔═╡ 9e6b9683-458c-4796-94b0-590b7325c11b
md"# Chapter 8: Matrix Norms and Perturbation Theory forLinear Systems"

# ╔═╡ 1ac0a8f6-6803-43bc-bc1e-8ac8d8b3ef1a
md"## 8.1 Vector Norms"

# ╔═╡ f613474b-8e87-4c3f-b352-570543e3990c
let
	x = [-2;5;4]
	my_norm(x,p)=sum(map(y->abs(y)^p,x))^(1/p)
	norm(x,1000)
end

# ╔═╡ 4a6c39bf-04b0-4758-ba67-8648ed21817a
md"##  8.2 Matrix Norms"

# ╔═╡ 17d84f46-a294-4c12-b706-72ebc15a98e0
cm"""
The follwoing are matrix norms
```math
\|\boldsymbol{A}\|_S:=\sum_{i=1}^m \sum_{i=1}^n\left|a_{i j}\right|,\quad \textbf{sum norm}
```
```math
\|\boldsymbol{A}\|_F:=\left(\sum_{i=1}^m \sum_{i=1}^n\left|a_{i j}\right|^2\right)^{1 / 2}, \quad \textbf{Frobenius norm}
```
```math
\|\boldsymbol{A}\|_M:=\max _{i, j}\left|a_{i j}\right|.\quad \textbf{max norm}
```
"""

# ╔═╡ d595178d-1edf-45d4-b97f-d8c98a8c88a4
let
	A = rand([-1.0,2,2,-2,0,3,3,3,3,4,-4,4],4,4)
	sum(abs.(A))== norm(A,1)
	
end

# ╔═╡ ba7cc1b5-e6bd-4b8b-adac-fa1636927d5e


# ╔═╡ e291e3fb-3301-44e9-8a54-842ccba1cf20
md"### Consistent and Subordinate Matrix Norms"

# ╔═╡ 183a420c-cac8-47d1-b9ae-dff435f10858
md"### Operator Norms"

# ╔═╡ a1de0642-7ef4-4d4b-ad3b-c2829cda3393
cm"""
- ``\|I\|_F=?`` so is it an operator norm?
"""

# ╔═╡ bef48346-8bd5-4f24-9a1e-cc91f52318d1
md"###  The Operator ``p``-Norms"

# ╔═╡ 1beb011b-26bf-4718-85b1-c16747e184cd
let
	A=(1//15)*[14  4  16 
	2  22  13]
	sum(map(abs,A))
	sum(map(x->x^2,A))
	maximum(sum(eachcol(A)))
end

# ╔═╡ aa8a8293-eeb8-4327-bfe7-c3d131a506fe
md"###  Unitary Invariant Matrix Norms"

# ╔═╡ b7bd9f6d-8df0-4cd4-b51c-2585da332d4c
let
	A = rand(ComplexF64,4,3)
	B = A'
	D = rand(ComplexF64,4,3)
	Q,R = qr(D)
	U=collect(Q)[1:3,:]
	V=collect(Q)[2:end,:]
	norm(V*A*U,2),norm(A,2),norm(collect(Q),2)
end

# ╔═╡ 5a1b36d2-3065-4fb2-85e3-c5f131f2cd26
md"### Absolute and Monotone Norms"

# ╔═╡ 7cb3af49-e0c4-4692-9efd-605d53189a54
md"## 8.3 The Condition Number with Respect to Inversion"

# ╔═╡ ef6598c7-2737-424a-8bbb-51298c9421a0
cm"""
- The number ``K(\boldsymbol{A})`` is called the __condition number with respect to inversion of a matrix__, or just the __condition number__ of ``\boldsymbol{A}``, if it is clear from the context that we are talking about inverting a matrix. 
- The condition number depends on the matrix ``\boldsymbol{A}`` and on the norm used. 
- If ``K(\boldsymbol{A})`` is large, ``\boldsymbol{A}`` is called __ill-conditioned__ (with respect to inversion). 
- If ``K(\boldsymbol{A})`` is small, ``\boldsymbol{A}`` is called __well-conditioned__ (with respect to inversion). 
- We always have ``K(\boldsymbol{A}) \geq 1``. 

For since ``\|\boldsymbol{x}\|=\|\boldsymbol{I} \boldsymbol{x}\| \leq\|\boldsymbol{I}\|\|\boldsymbol{x}\|`` for any ``\boldsymbol{x}`` we have ``\|\boldsymbol{I}\| \geq 1`` and therefore ``\|\boldsymbol{A}\|\left\|\boldsymbol{A}^{-1}\right\| \geq`` ``\left\|\boldsymbol{A} \boldsymbol{A}^{-1}\right\|=\|\boldsymbol{I}\| \geq 1``.
"""

# ╔═╡ 3a32d50b-75c5-407d-b0c4-2c4866ce7797
md"# Chapter 9: Least Squares"

# ╔═╡ b3a6b457-b7ff-481f-9192-f421ba34dbf5
cm"""
The __normal equations__ are
```math
\begin{aligned}
& \boldsymbol{A}^* \boldsymbol{A} \boldsymbol{x}=\left[\begin{array}{ccc}
1 & \cdots & 1 \\
t_1 & \cdots & t_m
\end{array}\right]\left[\begin{array}{ll}
1 & t_1 \\
\vdots \\
1 & t_m
\end{array}\right]\left[\begin{array}{l}
x_1 \\
x_2
\end{array}\right]=\left[\begin{array}{cc}
m & \sum t_k \\
\sum t_k & \sum t_k^2
\end{array}\right]\left[\begin{array}{l}
x_1 \\
x_2
\end{array}\right], \\
& =\left[\begin{array}{ccc}
1 & \cdots & 1 \\
t_1 & \cdots & t_m
\end{array}\right]\left[\begin{array}{c}
y_1 \\
\vdots \\
y_m
\end{array}\right]=\left[\begin{array}{c}
\sum y_k \\
\sum t_k y_k
\end{array}\right]=\boldsymbol{A}^* \boldsymbol{b},
\end{aligned}
```
"""

# ╔═╡ d21788c6-14d2-4237-b09e-0d869819aaf1
let
	t = 1:1.0:4
	y =[3.1;1.8;1.0;0.1]
	p = scatter(t,y,frame_style=:origin;m=:star,xlims=(-2,6),ylims=(-1,6),annotations=[(6,0.3,L"t"),(0.3,5.7,L"y")])

	A =hcat(ones(4),collect(t))
	b = A'*y
	B=A'*A
	x = B\b
	plot(p,t->x[1]+x[2]*t)
end

# ╔═╡ 249f483c-30ba-4ed3-9497-f774646999fc
md"## 9.3 Numerical Solution"

# ╔═╡ 4ec999bf-a88f-40c8-92d6-0d578084ab5d
cm"""
Assume that ``m \geq n, \boldsymbol{A} \in \mathbb{C}^{m \times n}, \boldsymbol{b} \in \mathbb{C}^m``

__Numerical Methods:__
1. Normal Equations
2. ``QR`` factorization
3. ``SVD`` factorization
"""

# ╔═╡ 67977543-3768-478f-b09f-a176a80f3f24
md"### Normal Equations"

# ╔═╡ f29b912a-7080-47f4-93cf-239a692bf337
cm"""
Assume ``\text{rank}(A) = n`` , then ``B = A^*A`` is positive definite.

- Solve using Cholesky factorization.

We calculate the inner product in ``A^*A`` in two ways
1. __inner product__: 
```math 
\left(\boldsymbol{A}^* \boldsymbol{A}\right)_{i, j}=\sum_{k=1}^m \bar{a}_{k, i} a_{k, j}, i, j=1, \ldots, n,\quad \left(\boldsymbol{A}^* \boldsymbol{b}\right)_i=\sum_{k=1}^m \bar{a}_{k, i} b_k, i=1, \ldots, n,
```   
2. __outer product__: 
```math
\boldsymbol{A}^* \boldsymbol{A}=\sum_{k=1}^m\left[\begin{array}{c}\bar{a}_{k, 1} \\ \vdots \\ \bar{a}_{k, n}\end{array}\right]\left[a_{k 1} \cdots a_{k n}\right],\quad \boldsymbol{A}^* \boldsymbol{b}=\sum_{k=1}^m\left[\begin{array}{c}\bar{a}_{k, 1} \\ \vdots \\ \bar{a}_{k, n}\end{array}\right] b_k.
```
- The __outer product__ form is suitable for large problems since it uses only one pass through the data importing one row of ``\boldsymbol{A}`` at a time from some separate storage.


"""

# ╔═╡ ceaf76d6-d5ee-4334-a3d8-88373f5f0f31
md"### QR Factorization"

# ╔═╡ fb1e4d9c-cae8-4a42-b218-00f157ba7b60
cm"""
Suppose ``\boldsymbol{A}=\boldsymbol{Q}_1 \boldsymbol{R}_1``, then
```math
\boldsymbol{A}^* \boldsymbol{A}=\boldsymbol{R}_1^* \boldsymbol{Q}_1^* \boldsymbol{Q}_1 \boldsymbol{R}_1=\boldsymbol{R}_1^* \boldsymbol{R}_1, \quad \boldsymbol{A}^* \boldsymbol{b}=\boldsymbol{R}_1^* \boldsymbol{Q}_1^* \boldsymbol{b} .
```
Since ``\boldsymbol{A}`` has rank ``n`` the matrix ``\boldsymbol{R}_1^*`` is nonsingular and can be canceled. Thus
```math
A^* A \boldsymbol{x}=\boldsymbol{A}^* \boldsymbol{b} \Longrightarrow \boldsymbol{R}_1 \boldsymbol{x}=\boldsymbol{c}_1, \quad \boldsymbol{c}_1:=\boldsymbol{Q}_1^* \boldsymbol{b} .
```
"""

# ╔═╡ 85794fff-8d0d-4ca3-bf94-b2aead8c9dd3
TableOfContents(title="📚 MATH557: Applied Linear Algebra", indent=true,depth=4)

# ╔═╡ 6449b443-e4e4-40d9-aefc-a98d0dd65cea
let
	x= [3;3;-1;-1] 
	y= [1;7;-4;2] 
	A = [ones(4) x y]
	b =ones(4)
	Q, R = qr(A)
	Q=collect(Q)
	Q1 = Q[:,1:3]
	# R1*x = Q1^*b
	x = R\Q1'*b
	# x = [1, 0, 0]
	# A\b
end

# ╔═╡ a0c5f1c5-4230-4630-8b4a-5515e7c49ebe
md"###  Singular Value Decomposition, Generalized Inverses and  Least Squares"

# ╔═╡ 4a35b5c1-086f-4103-b8db-98886ddbbf9d
cm"""
The matrix ``\boldsymbol{A}^{\dagger}`` is called the generalized inverse of ``\boldsymbol{A}``. We note that
1. If ``\boldsymbol{A}`` is square and nonsingular then ``\boldsymbol{A}^{-1}`` satisfies (👽) so that ``\boldsymbol{A}^{-1}=\boldsymbol{A}^{\dagger}``. Indeed, ``\boldsymbol{A}^{-1} \boldsymbol{A}=\boldsymbol{A} \boldsymbol{A}^{-1}=\boldsymbol{I}`` implies that ``\boldsymbol{A}^{-1} \boldsymbol{A}`` and ``\boldsymbol{A} \boldsymbol{A}^{-1}`` are Hermitian. Moreover, ``\boldsymbol{A} \boldsymbol{A}^{-1} \boldsymbol{A}=\boldsymbol{A}, \boldsymbol{A}^{-1} \boldsymbol{A} \boldsymbol{A}^{-1}=\boldsymbol{A}^{-1}``. By uniqueness ``\boldsymbol{A}^{-1}=\boldsymbol{A}^{\dagger}``.
2. We show in Exercise 9.7 that if ``\boldsymbol{A}`` has linearly independent columns then
```math
\boldsymbol{A}^{\dagger}=\left(\boldsymbol{A}^* \boldsymbol{A}\right)^{-1} \boldsymbol{A}^*
```
"""

# ╔═╡ b84719e4-902e-4c13-a146-74c19c28ca5e
let
	x= [3;3;-1;-1] 
	y= [1;7;-4;2] 
	A = [ones(4) x y]
	b=ones(4)
	# x = A\b # x = [1;0;0]
	U,s,V = svd(A,full=true)
	S = [Diagonal(s);zeros(1,3)]
	S1 =S[1:3,:]
	U1=U[:,1:3]
	V1=V
	S1inv=Diagonal(1 ./s)
	Ap=V1*S1inv*U1'
	x = Ap*b
end

# ╔═╡ e61790b0-ca31-44f8-bc2b-2bf610966400
md"# Chapter 12: The Classical Iterative Methods"

# ╔═╡ 9c266479-0995-4f72-9942-0941bb0938fd
md"""
Suppose $\boldsymbol{A} \in \mathbb{C}^{n \times n}$ is nonsingular and $\boldsymbol{b} \in \mathbb{C}^n$.
"""

# ╔═╡ 1557cf00-e2db-49df-b1b7-07b83ff6e214
md"## Component Form"

# ╔═╡ 0edaff9d-8630-435d-b2c5-4afea2ad7e5f
md"""

 

Let $\boldsymbol{x}_k=\left[\boldsymbol{x}_k(1), \ldots, \boldsymbol{x}_k(n)\right]^T$ be an approximation to the exact solution $\boldsymbol{x}$ of $\boldsymbol{A x}=\boldsymbol{b}$. 

We need to assume that the rows are ordered so that $\boldsymbol{A}$ has nonzero

1. __Jacobi's method ( $\mathbf{J}$ method)__ 


$\boldsymbol{x}_{k+1}(i)=\left(-\sum_{j=1}^{i-1} a_{i j} \boldsymbol{x}_k(j)-\sum_{j=i+1}^n a_{i j} \boldsymbol{x}_k(j)+b_i\right) / a_{i i}, \text { for } i=1,2, \ldots, n$

2. __Gauss-Seidel's method (GS method)__ is a modification of Jacobi's method, where we use the new $\boldsymbol{x}_{k+1}(i)$ immediately after it has been computed.

$\boldsymbol{x}_{k+1}(i)=\left(-\sum_{j=1}^{i-1} a_{i j} \boldsymbol{x}_{k+1}(j)-\sum_{j=i+1}^n a_{i j} \boldsymbol{x}_k(j)+b_i\right) / a_{i i}, \text { for } i=1,2, \ldots, n$

3. __The Successive overrelaxation method (SOR method)__ is obtained by introducing an acceleration parameter $0<\omega<2$ in the GS method. We write $\boldsymbol{x}(i)=\omega \boldsymbol{x}(i)+(1-\omega) \boldsymbol{x}(i)$ and this leads to the method
$\boldsymbol{x}_{k+1}(i)=\omega\left(-\sum_{j=1}^{i-1} a_{i j} \boldsymbol{x}_{k+1}(j)-\sum_{j=i+1}^n a_{i j} \boldsymbol{x}_k(j)+b_i\right) / a_{i i}+(1-\omega) \boldsymbol{x}_k(i)$

we can write this as

$\boldsymbol{x}_{k+1}=\omega \boldsymbol{x}^{gs}_{k+1} +(1-\omega) \boldsymbol{x}_k$

4. __The symmetric successive overrelaxation method SSOR__: Two SOR sweeps

- __Forward Sweep__
$\boldsymbol{x}_{k+1/2}=\omega \boldsymbol{x}^{gs}_{k+1} +(1-\omega) \boldsymbol{x}_k$

- __Backward Sweep__ For $i=n, n-1, ...,1$
$\small\boldsymbol{x}_{k+1}(i)=\omega\left(-\sum_{j=1}^{i-1} a_{i j} \boldsymbol{x}_{k+1/2}(j)-\sum_{j=i+1}^n a_{i j} \boldsymbol{x}_{k+1}(j)+b_i\right) / a_{i i}+(1-\omega) \boldsymbol{x}_{k+1/2}(i)$
"""

# ╔═╡ 3134e281-e748-4637-b182-97d27a14955f
md"## Matrix-Form"

# ╔═╡ d09b1a7b-90fd-46aa-88b9-954dd01d8043
cm"""
We write ``\boldsymbol{A x}=\boldsymbol{b}`` in the equivalent form
```math
M x=(M-A) x+b .
```
where 
- ``\boldsymbol{M}`` is nonsingular matrix 
- The matrix ``\boldsymbol{M}`` is known as a splitting matrix.

The corresponding iterative method is given by
```math
\boldsymbol{M} \boldsymbol{x}_{k+1}=(\boldsymbol{M}-\boldsymbol{A}) \boldsymbol{x}_k+\boldsymbol{b}
```
or
```math
\boldsymbol{x}_{k+1}:=\boldsymbol{G} \boldsymbol{x}_k+\boldsymbol{c}, \quad \boldsymbol{G}=\boldsymbol{I}-\boldsymbol{M}^{-1} \boldsymbol{A}, \quad, \boldsymbol{c}=\boldsymbol{M}^{-1} b .
```

- This is known as a __fixed-point iteration__. 
- The matrix ``\boldsymbol{M}`` can also be interpreted as a preconditioning matrix and is chosen to reduce the condition number of ``\boldsymbol{A}``.
"""

# ╔═╡ 1de083fd-d518-4c98-b27a-54ff74b6eb25
md"### The Splitting Matrices for the Classical Methods"

# ╔═╡ a387c180-f6ef-4399-8df3-35c72e2330e2
cm"""
```math
A=D-A_L-A_R,
```
where 
```math
\boldsymbol{A}_{\boldsymbol{L}}:=\left[\begin{array}{ccc}0 & & \\ -a_{2,1} & 0 & \\ \vdots & \ddots & \ddots \\ -a_{n, 1} & \cdots & -a_{n, n-1}\end{array}\right], \quad \boldsymbol{A}_{\boldsymbol{R}}:=\left[\begin{array}{ccc}0-a_{1,2} & \cdots & -a_{1, n} \\ \ddots & \ddots & \vdots \\ & 0 & -a_{n-1, n} \\ & & 0\end{array}\right].
```
and 
```math
\boldsymbol{D}:=\operatorname{diag}\left(a_{11}, \ldots, a_{n n}\right)
```

"""

# ╔═╡ 3932efd6-8805-4b31-a724-0c45685fdb6f
# J method
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	b=[6;25;-11;15]
	n = length(b)
	x0=zeros(n)
	xk=zeros(n)
	
	for k=1:10
		for i in 1:n
			xk[i] =(sum(-A[i,j]*x0[j] for j in 1:n if i!=j)+b[i])/A[i,i]
		end
		x0=copy(xk)
	end
	xk
end

# ╔═╡ 315db979-0472-45eb-a20a-9964821809f3
# J method
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	b=[6;25;-11;15]
	n = length(b)
	x0=zeros(n)
	xk=zeros(n)
	D = Diagonal(diag(A))
	Al=-tril(A,-1)
	Au=-triu(A,1)
	Minv = inv(D)
	G =(I-Minv*A)
	c =Minv*b
	for k=1:10
		xk=G*x0 + c
		x0=copy(xk)
	end
	xk
	# E = eigen(G)
	# maximum(abs.(E.values))
end

# ╔═╡ 279aea14-a829-4528-a7f6-d84eda6891ed
# GS method
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	b=[6;25;-11;15]
	n = length(b)
	x0=zeros(n)
	xk=zeros(n)
	
	for k=1:5
		for i in 1:n
			s1 = if i==1 
				0 
			else 
				sum(-A[i,j]*xk[j] for j in 1:n if j<i)
			end
			s2 = if i==n 
				0 
			else 
				sum(-A[i,j]*x0[j] for j in 1:n if j>i)
			end
			xk[i] =(s1+s2+b[i])/A[i,i]
		end
		x0=copy(xk)
	end
	xk
end

# ╔═╡ 87a8de33-2a99-4340-ad36-3d39abcc4de3
# GS method
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	b=[6;25;-11;15]
	n = length(b)
	x0=zeros(n)
	xk=zeros(n)
	D = Diagonal(diag(A))
	Al=-tril(A,-1)
	Au=-triu(A,1)
	M = D-Al
	Minv = inv(M)
	G = (I-Minv*A)
	c = Minv*b
	for k=1:10
		xk=G*x0 + c 
		x0=copy(xk)
	end
	xk
	# E = eigen(G)
	# maximum(abs.(E.values))
end

# ╔═╡ 14f1e18c-c521-4084-a13e-74a5e54fe572
# SOR method
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	# A =[4 3 0;3 4 -1;0 -1 4]
	b=[6;25;-11;15]
	# b=[24;30;-24]
	n = length(b)
	x0=ones(n)
	xk=ones(n)
	ω = 1.25
	for k=1:3
		for i in 1:n
			s1 = if i==1 
				0 
			else 
				sum(-A[i,j]*xk[j] for j in 1:n if j<i)
			end
			s2 = if i==n 
				0 
			else 
				sum(-A[i,j]*x0[j] for j in 1:n if j>i)
			end
			xk[i] =ω*(s1+s2+b[i])/A[i,i] + (1-ω)*x0[i]
		end
		
		x0=copy(xk)
	end
	xk
end


# ╔═╡ ba55ba2a-cb4e-4b31-b0df-9cf6593d56bc
# SOR method
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	b=[6;25;-11;15]
	ω=1.25
	n = length(b)
	x0=zeros(n)
	xk=zeros(n)
	D = Diagonal(diag(A))
	Al=-tril(A,-1)
	Au=-triu(A,1)
	
	M = (1/ω)*D-Al
	Minv = inv(M)
	G = (I-Minv*A)
	c = Minv*b
	for k=1:4
		xk=G*x0 + c
		x0=copy(xk)
	end
	xk
	# norm(G)
	E = eigen(G)
	maximum(abs.(E.values))
end

# ╔═╡ abf0c10e-b7d2-4a19-bfd5-2625015bccef
let
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	# A =[4 3 0;3 4 -1;0 -1 4]
	b=[6;25;-11;15]
	# b=[24;30;-24]
	n = length(b)
	x0=ones(n)
	xk=ones(n)
	ω = 1.25
	for k=1:3
		for i in 1:n
			s1 = if i==1 
				0 
			else 
				sum(-A[i,j]*xk[j] for j in 1:n if j<i)
			end
			s2 = if i==n 
				0 
			else 
				sum(-A[i,j]*x0[j] for j in 1:n if j>i)
			end
			xk[i] =(s1+s2+b[i])/A[i,i]
		end
		xk2=copy(xk)
		for i in n:-1:1
			xk[i]=ω*xk2[i] + (1-ω)*x0[i]
		end
		x0=copy(xk)
	end
	xk
end

# ╔═╡ ddb1b45e-df33-445c-b4b9-7611bc6f88a2
md"## Convergence "

# ╔═╡ 9222bbd2-a65e-425c-ac45-2abad5f43a5a
md"##  Stopping Criteria"

# ╔═╡ a5d9fa3d-7715-4bc3-88ee-44871493c84d
cm"""
1. ``\|x_{k}-x_{k-1}\|``
2. ``\|G\|``
"""

# ╔═╡ 08b43888-cc2b-4b0c-88fd-7af121088db1
# J method
let
	function Jmethod(A,b,ϵ, maxiters=100)
	n = length(b)
	x0=zeros(n)
	xk=zeros(n)
	D = Diagonal(diag(A))
	Al=-tril(A,-1)
	Au=-triu(A,1)
	Minv = inv(D)
	G =(I-Minv*A)
	c =Minv*b
	k = 1
	norm_of_G = norm(G)
	while true
		if k>maxiters 
			return xk, k, "max itrs reached"
		end
		
		xk=G*x0 + c
		
		# if norm_of_G^(k) <= ϵ 
		# 	return xk,k, "converged G"
		# end
		if norm(xk-x0) <= ϵ 
			return xk,k, "converged"
		end
		
		x0=copy(xk)
		k = k + 1
	end
	
	end
	A =[10 -1 2 0;-1 11 -1 3;2 -1 10 -1;0 3 -1 8]
	b=[6;25;-11;15]
	x,k, message = Jmethod(A,b,1e-5)
	# E = eigen(G)
	# maximum(abs.(E.values))
end

# ╔═╡ ed7ac1ae-3da3-4a46-a34b-4b445d52a95f
initialize_eqref()

# ╔═╡ 7b9ffd7c-3b93-4cfd-bed5-1590368ce987

@htl("""
<style>
@import url("https://mmogib.github.io/math102/custom.css");

ul {
  list-style: none;
}

ul li:before {
  content: '💡 ';
}

.p40 {
	padding-left: 40px;
}
</style>
""")

# ╔═╡ d0060e13-0aa0-495e-828e-084df48ef7e7

begin
    struct LocalImage
        filename
    end

    function Base.show(io::IO, ::MIME"image/png", w::LocalImage)
        write(io, read(w.filename))
    end
end

# ╔═╡ a8e556ea-48ca-4a43-914a-7171a6491186
LocalImage("docs/imgs/mathmatize_qrcode.png")

# exportqrcode("https://www.mathmatize.com/c/1263/","docs/imgs/mathmatize_qrcode.png",version=1,width=8,pixels=500)

# ╔═╡ d779340e-4dab-45c1-b8df-c0bcbae32a90

begin
	function add_space(n=1)
		repeat("&nbsp;",n)
	end
    function post_img(img::String, w=500)
        res = Resource(img, :width => w)
        cm"""
      <div class="img-container">

      $(res)

      </div>"""
    end
    function poolcode()
        cm"""
      <div class="img-container">

      $(Resource("https://www.dropbox.com/s/cat9ots4ausfzyc/qrcode_itempool.com_kfupm.png?raw=1",:width=>300))

      </div>"""
    end
    function define(t="")
        beginBlock("Definition", t)
    end
    function bbl(t)
        beginBlock(t, "")
    end
    function bbl(t, s)
        beginBlock(t, s)
    end
    ebl() = endBlock()
	function theorem(s)
		bth(s)
	end
    function bth(s)
        beginTheorem(s)
    end
    eth() = endTheorem()
	ex() = example("Example","")
	ex(n::Int; s::String="") = ex("Example $n", s)
    ex(t, s) = example(t, s)
    function beginBlock(title, subtitle)
        """<div style="box-sizing: border-box;">
       	<div style="display: flex;flex-direction: column;border: 6px solid rgba(200,200,200,0.5);box-sizing: border-box;">
       	<div style="display: flex;">
       	<div style="background-color: #FF9733;
       	    border-left: 10px solid #df7300;
       	    padding: 5px 10px;
       	    color: #fff!important;
       	    clear: left;
       	    margin-left: 0;font-size: 112%;
       	    line-height: 1.3;
       	    font-weight: 600;">$title</div>  <div style="olor: #000!important;
       	    margin: 0 0 20px 25px;
       	    float: none;
       	    clear: none;
       	    padding: 5px 0 0 0;
       	    margin: 0 0 0 20px;
       	    background-color: transparent;
       	    border: 0;
       	    overflow: hidden;
       	    min-width: 100px;font-weight: 600;
       	    line-height: 1.5;">$subtitle</div>
       	</div>
       	<p style="padding:5px;">
       """
    end
    function beginTheorem(subtitle)
        beginBlock("Theorem", subtitle)
    end
    function endBlock()
        """</p></div></div>"""
    end
    function endTheorem()
        endBlock()
    end
    function example(lable, desc)
        """<div style="display:flex;">
       <div style="
       font-size: 112%;
           line-height: 1.3;
           font-weight: 600;
           color: #f9ce4e;
           float: left;
           background-color: #5c5c5c;
           border-left: 10px solid #474546;
           padding: 5px 10px;
           margin: 0 12px 20px 0;
           border-radius: 0;
       ">$lable:</div>
       <div style="flex-grow:3;
       line-height: 1.3;
           font-weight: 600;
           float: left;
           padding: 5px 10px;
           margin: 0 12px 20px 0;
           border-radius: 0;
       ">$desc</div>
       </div>"""
    end
    @htl("")
end

# ╔═╡ 06790633-b8aa-45ba-ac67-c416c88b166a
begin
    text_book = post_img("https://www.dropbox.com/scl/fi/0axh6gwcwp4eg97wxfwtd/text_book.webp?rlkey=zlptduxz777cvreoghlcerwbm&raw=1", 200)
    md""" # Syllabus
    ## Syallbus
    See here [Term 241 - MATH557 - Syllabus](https://www.dropbox.com/scl/fi/y7fnigu9r4ez5ilma11qq/T241_MATH557_Syllabus.pdf?rlkey=8w3lyyqvsn59eyp1s0pwcslbk&raw=1)
    ## Textbook
    __Textbook: Lyche, T. (2020). Numerical Linear Algebra and Matrix Factorizations (Vol. 22). Springer International Publishing. [https://doi.org/10.1007/978-3-030-36468-7](https://doi.org/10.1007/978-3-030-36468-7)__
    $text_book

    ## Office Hours
    I strongly encourage all students to make use of my office hours. These dedicated times are a valuable opportunity for you to ask questions, seek clarification on lecture material, discuss challenging problems, and get personalized feedback on your work. Engaging with me during office hours can greatly enhance your understanding of the course content and improve your performance. Whether you're struggling with a specific concept or simply want to delve deeper into the subject, I am here to support your learning journey. Don't hesitate to drop by; __your success is my priority__.

    | Day       | Time        |
    |-----------|-------------|
    | Tuesday    | 08:00-08:40PM |
    | Thursday | 08:00-08:40PM |
    Also you can ask for an online meeting through __TEAMS__.
    """
end

# ╔═╡ c6241f8a-60ff-44c0-b363-f2d91b2c5eb0
cm"""
__Applied/Numerical linear algebra__ providesalgorithms for solving problems of the
following kind:

$(bbl("System of linear equations",""))

Given a (square) matrix ``\boldsymbol{A}`` and a vector ``\boldsymbol{b}``. 

Find a vector ``\boldsymbol{x}`` such that ``\boldsymbol{A x}=\boldsymbol{b}``.
$(ebl())

$(bbl("Least squares",""))

Given a (rectangular) matrix ``\boldsymbol{A}`` and a vector ``\boldsymbol{b}``. 

Find a vector ``\boldsymbol{x}`` such that the sum of squares of the components of ``\boldsymbol{b}-\boldsymbol{A} \boldsymbol{x}`` is as small as possible.
$(ebl())

$(bbl("Eigenvalues and eigenvectors",""))
Given a (square) matrix ``\boldsymbol{A}``. 

Find a number ``\lambda`` and/or a nonzero vector ``\boldsymbol{x}`` such that ``\boldsymbol{A} \boldsymbol{x}=\lambda \boldsymbol{x}``.
$(ebl())
"""

# ╔═╡ 6dec69df-0c80-4790-b5c5-12fbe2dc41b8
cm"""
Associated with an ``m \times n`` matrix ``\boldsymbol{X}=\left[\boldsymbol{x}_1, \ldots, \boldsymbol{x}_n\right]``, where ``\boldsymbol{x}_j \in \mathcal{V}, j=1, \ldots, n`` are the following subspaces of ``\mathcal{V}``.
1. The subspace 
```math
\mathcal{R}(\boldsymbol{X}):=\left\{\boldsymbol{X} \boldsymbol{c}: \boldsymbol{c} \in \mathbb{R}^n\right\}=\operatorname{span}(\mathcal{X})
``` 
$(add_space(10))is called the column space of ``\boldsymbol{X}``. It is the smallest subspace containing ``\mathcal{X}=\left\{\boldsymbol{x}_1, \ldots, \boldsymbol{x}_n\right\}``. 
$(add_space(20))▶ The dimension of ``\mathcal{R}(\boldsymbol{X})`` is called the __rank of ``\boldsymbol{X}``__. 

$(add_space(20))▶ The matrix ``\boldsymbol{X}`` has rank ``n`` if and only if it has linearly independent columns.
{}
2. ``\mathcal{R}\left(\boldsymbol{X}^T\right)`` is called the row space of ``\boldsymbol{X}``. It is generated by the rows of ``\boldsymbol{X}`` written as column vectors.
3. The subspace ``\mathcal{N}(\boldsymbol{X}):=\left\{\boldsymbol{y} \in \mathbb{R}^n: \boldsymbol{X} \boldsymbol{y}=\boldsymbol{0}\right\}`` is called the __null space__ or __kernel space__ of ``\boldsymbol{X}``. The dimension of ``\mathcal{N}(\boldsymbol{X})`` is called the __nullity__ of ``\boldsymbol{X}`` and denoted ``\operatorname{null}(\boldsymbol{X})``.

4. Suppose ``\boldsymbol{X} \in`` ``\mathbb{C}^{m \times n}``. Then
	- ``\operatorname{rank}(X)=\operatorname{rank}\left(X^*\right)`` where ``X^*=\overline{X}^T``.
	- ``\operatorname{rank}(\boldsymbol{X})+\operatorname{null}(\boldsymbol{X})=n``,
	- ``\operatorname{rank}(\boldsymbol{X})+\operatorname{null}\left(\boldsymbol{X}^*\right)=m``,
"""

# ╔═╡ 5742ff61-39a1-41e6-9136-1f148516d5e0
cm"""
The system is
- __homogeneous__ if ``\boldsymbol{b}=\boldsymbol{0}`` and 
- __underdetermined__ $(add_space(40)) if $(add_space(10)) ``m< n``
- __square__ $(add_space(67)) if $(add_space(10)) ``m=n``
- __overdetermined__ $(add_space(44)) if $(add_space(10))``m > n``

The matrix ``\boldsymbol{A} \in\mathbb{R}^{n \times n}`` is called
- __nonsingular__ if the only real solution of the homogeneous system ``\boldsymbol{A} \boldsymbol{x}=\mathbf{0}`` is ``\boldsymbol{x}=\mathbf{0}``. 
- __singular__ if there is a nonzero ``\boldsymbol{x} \in \mathbb{R}^n`` such that ``\boldsymbol{A x}=\mathbf{0}``.
"""

# ╔═╡ 7f18dc1e-e040-4f0d-b3bf-a5477489a1ab
cm"""
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is a square matrix. 

- A matrix ``\boldsymbol{B} \in \mathbb{C}^{n \times n}`` is called a __right inverse of__ ``\boldsymbol{A}`` if 
```math
\boldsymbol{A B}=\boldsymbol{I}.
``` 
- A matrix ``\boldsymbol{C} \in \mathbb{C}^{n \times n}`` is said to be a __left inverse__ of ``\boldsymbol{A}`` if ``\boldsymbol{C} \boldsymbol{A}=\boldsymbol{I}``. 
- We say that ``\boldsymbol{A}`` is __invertible__ if it has both a left and a right inverse. If ``\boldsymbol{A}`` has a right inverse ``\boldsymbol{B}`` and a left inverse ``\boldsymbol{C}`` then
```math
\boldsymbol{C}=\boldsymbol{C} \boldsymbol{I}=\boldsymbol{C}(\boldsymbol{A B})=(\boldsymbol{C} \boldsymbol{A}) \boldsymbol{B}=\boldsymbol{I} \boldsymbol{B}=\boldsymbol{B}.
```
$(add_space(10))It is called the inverse of ``\boldsymbol{A}`` and denoted by ``\boldsymbol{A}^{-1}``. 
- Thus the inverse satisfies ``\boldsymbol{A}^{-1} \boldsymbol{A}=\boldsymbol{A} \boldsymbol{A}^{-1}=\boldsymbol{I}``.
-  A square matrix is invertible if and only if it is nonsingular.
Suppose ``\boldsymbol{A}, \boldsymbol{B} \in \mathbb{C}^{n \times n}`` are nonsingular and ``c`` is a nonzero constant.
- ``\boldsymbol{A}^{-1}`` is nonsingular and ``\left(\boldsymbol{A}^{-1}\right)^{-1}=\boldsymbol{A}``.
- ``\boldsymbol{C}=\boldsymbol{A} \boldsymbol{B}`` is nonsingular and ``\boldsymbol{C}^{-1}=\boldsymbol{B}^{-1} \boldsymbol{A}^{-1}``.
- ``\boldsymbol{A}^T`` is nonsingular and ``\left(\boldsymbol{A}^T\right)^{-1}=\left(\boldsymbol{A}^{-1}\right)^T=: \boldsymbol{A}^{-T}``.
- ``\boldsymbol{A}^*`` is nonsingular and ``\left(\boldsymbol{A}^*\right)^{-1}=\left(\boldsymbol{A}^{-1}\right)^*=: \boldsymbol{A}^{-*}``.
- ``c \boldsymbol{A}`` is nonsingular and ``(c \boldsymbol{A})^{-1}=\frac{1}{c} \boldsymbol{A}^{-1}``.
"""

# ╔═╡ 0fb73bc3-e70d-4de2-a8fa-49938c6f3a60
cm"""
$(bbl("Properties",""))

- if ``A`` is triangular then ``\operatorname{det}(A) = a_{11} a_{22} \cdots a_{n n}``, the product of the diagonal elements. 
- ``\operatorname{det}(\boldsymbol{I})=1``. 

The elementary operations using either rows or columns are
- Interchanging two rows(columns): ``\operatorname{det}(\boldsymbol{B})=-\operatorname{det}(\boldsymbol{A})``,
- Multiply a row(column) by a scalar: ``\alpha, \operatorname{det}(\boldsymbol{B})=\alpha \operatorname{det}(\boldsymbol{A})``,
- Add a constant multiple of one row(column) to another row(column): ``\operatorname{det}(\boldsymbol{B})=\operatorname{det}(\boldsymbol{A})``.

where ``\boldsymbol{B}`` is the result of performing the indicated operation on ``\boldsymbol{A}``.

```math
\begin{aligned}
& \operatorname{det}(\boldsymbol{A})=\sum_{j=1}^n(-1)^{i+j} a_{i j} \operatorname{det}\left(\boldsymbol{A}_{i j}\right) \text { for } i=1, \ldots, n, \text { row } \\
\end{aligned}
```

Here ``\boldsymbol{A}_{i, j}`` denotes the submatrix of ``\boldsymbol{A}`` obtained by deleting the ``i`` th row and ``j`` th column of ``\boldsymbol{A}``. 

For ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` and ``1 \leq i, j \leq n`` the determinant ``\operatorname{det}\left(\boldsymbol{A}_{i j}\right)`` is called the __cofactor__ of ``a_{i j}``.

-----------

If ``\boldsymbol{A}, \boldsymbol{B}`` are square matrices of order ``n`` with real or complex elements, then
- ``\operatorname{det}(\boldsymbol{A B})=\operatorname{det}(\boldsymbol{A}) \operatorname{det}(\boldsymbol{B})``.
- ``\operatorname{det}\left(\boldsymbol{A}^T\right)=\operatorname{det}(\boldsymbol{A})``, and ``\operatorname{det}\left(\boldsymbol{A}^*\right)=\overline{\operatorname{det}(\boldsymbol{A})}``, (complex conjugate).
- ``\operatorname{det}(a \boldsymbol{A})=a^n \operatorname{det}(\boldsymbol{A})``, for ``a \in \mathbb{C}``.
- ``\boldsymbol{A}`` is singular if and only if ``\operatorname{det}(\boldsymbol{A})=0``.
- If ``\boldsymbol{A}=\left[\begin{array}{ll}\boldsymbol{C} & \boldsymbol{D} \\ \mathbf{0} & \boldsymbol{E}\end{array}\right]`` for some square matrices ``\boldsymbol{C}, \boldsymbol{E}`` then ``\operatorname{det}(\boldsymbol{A})=\operatorname{det}(\boldsymbol{C}) \operatorname{det}(\boldsymbol{E})``.
- __Cramer's rule__ Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is nonsingular and ``\boldsymbol{b} \in \mathbb{C}^n``. Let ``\boldsymbol{x}=`` ``\left[x_1, x_2, \ldots, x_n\right]^T`` be the unique solution of ``\boldsymbol{A} \boldsymbol{x}=\boldsymbol{b}``. Then
```math
x_j=\frac{\operatorname{det}\left(\boldsymbol{A}_j(\boldsymbol{b})\right)}{\operatorname{det}(\boldsymbol{A})}, \quad j=1,2, \ldots, n
```
$(add_space(10))where ``\boldsymbol{A}_j(\boldsymbol{b})`` denote the matrix obtained from ``\boldsymbol{A}`` by replacing the ``j`` th column of ``\boldsymbol{A}`` by ``\boldsymbol{b}``.
- __Adjoint formula for the inverse.__ If ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is nonsingular then
```math
\boldsymbol{A}^{-1}=\frac{1}{\operatorname{det}(\boldsymbol{A})} \operatorname{adj}(\boldsymbol{A})
```

<div style="padding-left:40px;">

where the matrix ``\operatorname{adj}(\boldsymbol{A}) \in \mathbb{C}^{n \times n}`` with elements ``\operatorname{adj}(\boldsymbol{A})_{i, j}=(-1)^{i+j} \operatorname{det}\left(\boldsymbol{A}_{j, i}\right)`` is called the adjoint of ``\boldsymbol{A}``. Moreover, ``\boldsymbol{A}_{j, i}`` denotes the submatrix of ``\boldsymbol{A}`` obtained by deleting the ``j`` th row and ``i`` th column of ``\boldsymbol{A}``.

</div>

"""

# ╔═╡ 0beeed3b-905e-4ac1-b19c-a06cec3f160f
cm"""
📑 __finite difference method__ can use to approximate the solution.
```math
\begin{array}{lcl}
u'(x) &=& \displaystyle\lim_{h\to 0} \frac{u(x+h/2)-u(x-h/2)}{h}\\
u''(x) &=&\displaystyle \lim_{h\to 0} \frac{u(x+h)-2u(x)+u(x-h)}{h^2}\\
\end{array}
```
Let ``m`` be a positive integer, and 
- ``h:=1 /(m+1)`` be the discretization parameter, and 
- replace the interval ``[0,1]`` by grid points 
```math
x_j:=j h \quad \text{  for  } \quad j=0,1, \ldots, m+1. 
```
We then obtain approximations ``v_j`` to the exact solution ``u\left(x_j\right)`` for ``j=1, \ldots, m`` by replacing the differential equation by the difference equation
```math
\frac{-v_{j-1}+2 v_j-v_{j+1}}{h^2}=f(j h), \quad j=1, \ldots, m, \quad v_0=v_{m+1}=0
```

Moving the ``h^2`` factor to the right hand side this can be written as an ``m \times m`` linear system
$(post_img("https://www.dropbox.com/scl/fi/d2wh1w0zrcuid54h3euh3/eq_2_21.png?rlkey=7zdri03ldvs6si9a1iyf4w8vb&dl=1"))

- The matrix ``\boldsymbol{T}`` is called the __second derivative matrix__ and will occur frequently in this course. 
```math
\boldsymbol{T}=\operatorname{tridiag}\left(a_i, d_i, c_i\right) \in \mathbb{R}^{m \times m},
```
$(add_space(10))where in this case ``a_i=c_i=-1`` and ``d_i=2`` for all ``i``.
"""

# ╔═╡ 7513851f-48f8-4e2f-8c97-baee10e9e4c7
cm"""
$(bbl("Definition","(Diagonal Dominance)"))
The matrix ``\boldsymbol{A}=\left[a_{i j}\right] \in \mathbb{C}^{n \times n}`` is __weakly diagonally dominant__ if
```math
\left|a_{i i}\right| \geq \sum_{j \neq i}\left|a_{i j}\right|, i=1, \ldots, n
```
$(ebl())
- If ``\boldsymbol{A}=\operatorname{tridiag}\left(a_i, d_i, c_i\right) \in`` ``\mathbb{C}^{n \times n}`` is tridiagonal and weakly diagonally dominant. If in addition ``\left|d_1\right|>\left|c_1\right|`` and ``a_i \neq 0`` for ``i=1, \ldots, n-2``, then ``\boldsymbol{A}`` has a unique ``L U`` factorization $(eqref("lufactor")). If in addition ``d_n \neq 0``, then ``\boldsymbol{A}`` is nonsingular.
"""

# ╔═╡ 3dcc9b0b-2daa-4695-9e77-311672ac511b
cm"""
##### Properties of Block Multiplication

Assume that ``\boldsymbol{A} \in \mathbb{C}^{m \times p}`` and ``\boldsymbol{B} \in \mathbb{C}^{p \times n}``. Then

1. If ``\boldsymbol{B}=\left[\boldsymbol{b}_{: 1}, \ldots, \boldsymbol{b}_{: n}\right]`` is partitioned into columns then the partition of the product ``\boldsymbol{A} \boldsymbol{B}`` into columns is
```math
\boldsymbol{A B}=\left[\boldsymbol{A} \boldsymbol{b}_{: 1}, \boldsymbol{A} \boldsymbol{b}_{: 2}, \ldots, \boldsymbol{A} \boldsymbol{b}_{: n}\right]
```
$(add_space(12))In particular, if ``\boldsymbol{I}`` is the identity matrix of order ``p`` then
```math
\boldsymbol{A}=\boldsymbol{A} \boldsymbol{I}=\boldsymbol{A}\left[\boldsymbol{e}_1, \boldsymbol{e}_2, \ldots, \boldsymbol{e}_p\right]=\left[\boldsymbol{A} \boldsymbol{e}_1, \boldsymbol{A} \boldsymbol{e}_2, \ldots, \boldsymbol{A} \boldsymbol{e}_p\right]
```
$(add_space(12))and we see that column ``j`` of ``\boldsymbol{A}`` can be written ``\boldsymbol{A} \boldsymbol{e}_j`` for ``j=1, \ldots, p``.

2. Similarly, if ``\boldsymbol{A}`` is partitioned into rows then
```math
\boldsymbol{A} \boldsymbol{B}=\left[\begin{array}{c}
a_{1:}^T \\
a_{2:}^T \\
\vdots \\
a_{m:}^T
\end{array}\right] \boldsymbol{B}=\left[\begin{array}{c}
a_{1:}^T \boldsymbol{B} \\
a_{2:}^T \boldsymbol{B} \\
\vdots \\
a_{m:}^T \boldsymbol{B}
\end{array}\right]
```

3. It is often useful to write the matrix-vector product ``\boldsymbol{A x}`` as a linear combination of the columns of ``\boldsymbol{A}``
```math
\boldsymbol{A x}=x_1 \boldsymbol{a}_{: 1}+x_2 \boldsymbol{a}_{: 2}+\cdots+x_p \boldsymbol{a}_{: p}
```
"""

# ╔═╡ 3f43de8a-b6b4-4d53-aad4-bd3356f2d2a8
cm"""
8. Consider finally the general case. If all the matrix products ``\boldsymbol{A}_{i k} \boldsymbol{B}_{k j}`` in
```math
\boldsymbol{C}_{i j}=\sum_{k=1}^s \boldsymbol{A}_{i k} \boldsymbol{B}_{k j}, \quad i=1, \ldots, p, j=1, \ldots, q
```
$(add_space(10))are well defined then
```math
\left[\begin{array}{cccc}
\boldsymbol{A}_{11} & \cdots & \boldsymbol{A}_{1 s} \\
\vdots & & \vdots \\
\boldsymbol{A}_{p 1} & \cdots & \boldsymbol{A}_{p s}
\end{array}\right]\left[\begin{array}{ccc}
\boldsymbol{B}_{11} & \cdots & \boldsymbol{B}_{1 q} \\
\vdots & & \vdots \\
\boldsymbol{B}_{s 1} & \cdots & \boldsymbol{B}_{s q}
\end{array}\right]=\left[\begin{array}{ccc}
\boldsymbol{C}_{11} & \cdots & \boldsymbol{C}_{1 q} \\
\vdots & & \vdots \\
\boldsymbol{C}_{p 1} & \cdots & \boldsymbol{C}_{p q}
\end{array}\right]
```

##### The requirements are that
- the number of columns in ``\boldsymbol{A}`` is equal to the number of rows in ``\boldsymbol{B}``.
- the position of the vertical partition lines in ``\boldsymbol{A}`` has to mach the position of the horizontal partition lines in ``\boldsymbol{B}``. The horizontal lines in ``\boldsymbol{A}`` and the vertical lines in ``\boldsymbol{B}`` can be anywhere.
"""

# ╔═╡ 29f4fc5f-9f93-4004-9089-061fc14517f9
cm"""
$(bbl("Inverse of a Block Triangular Matrix")) Suppose
```math
\boldsymbol{A}=\left[\begin{array}{cc}
\boldsymbol{A}_{11} & \boldsymbol{A}_{12} \\
\mathbf{0} & \boldsymbol{A}_{22}
\end{array}\right]
```
where ``\boldsymbol{A}, \boldsymbol{A}_{11}`` and ``\boldsymbol{A}_{22}`` are square matrices. Then ``\boldsymbol{A}`` is nonsingular if and only if both ``\boldsymbol{A}_{11}`` and ``\boldsymbol{A}_{22}`` are nonsingular. In that case
```math
\boldsymbol{A}^{-1}=\left[\begin{array}{cc}
\boldsymbol{A}_{11}^{-1} & \boldsymbol{C} \\
\mathbf{0} & \boldsymbol{A}_{22}^{-1}
\end{array}\right]
```
for some matrix ``\boldsymbol{C}``.
"""

# ╔═╡ 1740ab31-1c55-4bfc-a34e-dc10852de7ff
cm"""
$(bbl("Inverse of a Triangular Matrix")) An upper (lower) triangular matrix ``\boldsymbol{A}=\left[a_{i j}\right] \in \mathbb{C}^{n \times n}`` is nonsingular if and only if the diagonal elements ``a_{i i}``, ``i=1, \ldots, n`` are nonzero. In that case the inverse is upper (lower) triangular with diagonal elements ``a_{i i}^{-1}, i=1, \ldots, n``.
"""

# ╔═╡ 45ef1c9d-1052-4656-b4e5-1704c39013ee
cm"""
$(bbl("Product of Triangular Matrices")) The product ``\boldsymbol{C}=\boldsymbol{A} \boldsymbol{B}=\left(c_{i j}\right)`` of two upper (lower) triangular matrices ``\boldsymbol{A}=\left(a_{i j}\right)`` and ``\boldsymbol{B}=\left(b_{i j}\right)`` is upper (lower) triangular with diagonal elements ``c_{i i}=a_{i i} b_{i i}`` for all ``i``.
"""

# ╔═╡ cda6aea6-61ca-483f-98f9-2274ef50e049
cm"""
$(bbl("Unit Triangular Matrices")) For a unit upper (lower) triangular (i.e. ``1'``s on the diagonal) matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` :
1. ``\boldsymbol{A}`` is nonsingular and the inverse is unit upper(lower) triangular.
2. The product of two unit upper (lower) triangular matrices is unit upper (lower) triangular.
"""

# ╔═╡ 58451143-9990-4426-9d42-50a5d3d9849b
cm"""
$(ex(1))
Using Gaussian Elemination (without row exchange) solve the system
```math
\begin{array}{lcllcllcllcllcllcllcl} 
x_1&+&x_2&&&+&3 x_4&= & 4 \\ 
2 x_1&+&x_2&-&x_3&+&x_4&=& 1 \\ 
3 x_1&-&x_2&-&x_3&+&2 x_4&=& -3 \\ 
-x_1&+&2 x_2&+&3 x_3&-&x_4&= & 4
\end{array}
```
"""

# ╔═╡ 5fef6c3a-9809-4906-b5b8-4d00c72439a8
cm"""
$(bbl("Remarks",""))
- Gauss elimination leads to an ``LU`` __factorization__. 
- Solving ``Ly=b`` is called __forward substitution__.
- Solving ``Uy=y`` is called __backward substitution__.
- 
"""

# ╔═╡ b5c0a76a-dc08-46f9-9b2d-f1b528f1aa87
cm"""
The process transforming ``\boldsymbol{A}^{(k)}`` into ``\boldsymbol{A}^{(k+1)}`` for ``k=1, \ldots, n-1`` can be described as follows.
```math
\begin{aligned}
& \text { for } i=k+1: n \\
& l_{i k}^{(k)}=a_{i k}^{(k)} / a_{k k}^{(k)} \\
& \text { for } j=k: n \\
& \quad a_{i j}^{(k+1)}=a_{i j}^{(k)}-l_{i k}^{(k)} a_{k j}^{(k)}
\end{aligned}
```

$(post_img("https://www.dropbox.com/scl/fi/p7pgghrjx7pn7xra4ks6k/fig3_1.png?rlkey=k12xy1j7quuf4tugjvwfcvl0k&dl=1",500))
"""

# ╔═╡ e141f527-3783-4293-a890-10c05aa166dd
cm"""
- ``\boldsymbol{A}^{(k+1)}`` will have zeros under the diagonal in its first ``k`` columns and the elimination is carried one step further. 
- The numbers ``l_{i k}^{(k)}`` in (3.2) are called multipliers.

__Gaussian elimination with no row interchanges is valid if and only if the pivots ``a_{k k}^{(k)}`` are nonzero for ``k=1, \ldots, n-1``__. 

$(define("3.1 (Principal Submatrix)")) 
For ``k=1, \ldots, n`` the matrices ``\boldsymbol{A}_{[k]} \in \mathbb{C}^{k \times k}`` given by
```math
\boldsymbol{A}_{[k]}:=\boldsymbol{A}(1: k, 1: k)=\left[\begin{array}{ccc}
a_{11} & \cdots & a_{k 1} \\
\vdots & & \vdots \\
a_{k 1} & \cdots & a_{k k}
\end{array}\right]
```
are called the __leading principal submatrices__ of ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``. More generally, a matrix ``\boldsymbol{B} \in \mathbb{C}^{k \times k}`` is called a __principal submatrix__ of ``\boldsymbol{A}`` if ``\boldsymbol{B}=\boldsymbol{A}(\boldsymbol{r}, \boldsymbol{r})``, where ``\boldsymbol{r}=\left[r_1, \ldots, r_k\right]`` for some ``1 \leq r_1< \cdots< r_k \leq n``. Thus,
```math
b_{i, j}=a_{r_i, r_j}, \quad i, j=1, \ldots, k
```

The determinant of a (leading) principal submatrix is called a __(leading) principal minor__.
"""

# ╔═╡ 6d5bf021-a49e-4726-a446-900422cc6703
cm"""
$(bth("3.1")) We have ``a_{k, k}^{(k)} \neq 0`` for ``k=1, \ldots, n-1`` if and only if the leading principal submatrices ``\boldsymbol{A}_{[k]}`` of ``\boldsymbol{A}`` are nonsingular for ``k=1, \ldots, n-1``. Moreover
```math
\operatorname{det}\left(\boldsymbol{A}_{[k]}\right)=a_{11}^{(1)} a_{22}^{(2)} \cdots a_{k k}^{(k)}, \quad k=1, \ldots, n
```
"""

# ╔═╡ a4226ef9-8d8f-4508-9e88-6dc09f6410be
cm"""
$(bth("3.2"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` and that the leading principal submatrices ``\boldsymbol{A}_{[k]}`` are nonsingular for ``k=1, \ldots, n-1``. Then Gaussian elimination with no row interchanges results in an ``L U`` factorization of ``\boldsymbol{A}``. In particular ``\boldsymbol{A}=\boldsymbol{L} \boldsymbol{U}``, where
```math
\boldsymbol{L}=\left[\begin{array}{cccc}
1 & & & \\
l_{21}^{(1)} & 1 & & \\
\vdots & & \ddots \\
l_{n 1}^{(1)} & l_{n 2}^{(2)} & \cdots & 1
\end{array}\right], \quad \boldsymbol{U}=\left[\begin{array}{ccc}
a_{11}^{(1)} & \cdots & a_{1 n}^{(1)} \\
& \ddots & \vdots \\
& & a_{n n}^{(n)}
\end{array}\right]
```
where the ``l_{i j}^{(j)}`` and ``a_{i j}^{(i)}`` are given by 
```math
\begin{aligned}
 l_{i j}^{(j)}&=a_{i j}^{(j)} / a_{j j}^{(j)} \\
 \quad a_{i j}^{(j+1)}&=a_{i j}^{(j)}-l_{i j}^{(j)} a_{j j}^{(j)}
\end{aligned}
```
"""

# ╔═╡ 6897cf15-2370-4114-a574-ea137096a175
cm"""
- Any nonsingular linear system can be solved by Gaussian elimination if we incorporate row interchanges
- Interchanging two rows(and/or two columns) during Gaussian elimination is known
 as __pivoting__. 

$(define("3.3")) A permutation matrix is a matrix of the form
```math
\boldsymbol{P}=\boldsymbol{I}(:, \boldsymbol{p})=\left[\boldsymbol{e}_{i_1}, \boldsymbol{e}_{i_2}, \ldots, \boldsymbol{e}_{i_n}\right] \in \mathbb{R}^{n \times n}
```
where ``\boldsymbol{e}_{i_1}, \ldots, \boldsymbol{e}_{i_n}`` is a permutation of the unit vectors ``\boldsymbol{e}_1, \ldots, \boldsymbol{e}_n \in \mathbb{R}^n``.
"""

# ╔═╡ a0560ea2-1ab5-425a-a418-92e196b51f92
cm"""
$(define(3.4))
We define a ``(j, k)``-Interchange matrix ``\boldsymbol{I}_{j k}`` by interchanging column ``j`` and ``k`` of the identity matrix.
"""

# ╔═╡ 82c571b9-23a6-48b6-b58c-a42269bfa429
cm"""
$(bth("3.3"))
Gaussian elimination with row pivoting on a nonsingular matrix ``\boldsymbol{A} \in`` ``\mathbb{C}^{n \times n}`` leads to the factorization ``\boldsymbol{A}=\boldsymbol{P} \boldsymbol{L} \boldsymbol{U}``, where ``\boldsymbol{P}`` is a permutation matrix, ``\boldsymbol{L}`` is lower triangular with ones on the diagonal, and ``\boldsymbol{U}`` is upper triangular. More explicitly, ``\boldsymbol{P}=\boldsymbol{I}(:, \boldsymbol{p})``, where ``\boldsymbol{p}=\boldsymbol{I}_{r_{n-1}, n-1} \cdots \boldsymbol{I}_{r_1, 1}[1, \ldots, n]^T``, and
```math
\boldsymbol{L}=\left[\begin{array}{ccc}
1 & & \\
a_{p_2, 1}^{(1)} & 1 & \\
\vdots & & \ddots \\
a_{p_n, 1}^{(1)} & a_{p_n, 2}^{(2)} & \cdots
\end{array}\right], \quad \boldsymbol{U}=\left[\begin{array}{ccc}
a_{p_1, 1}^{(1)} & \cdots & a_{p_1, n}^{(1)} \\
& \ddots & \vdots \\
& & a_{p_n, n}^{(n)}
\end{array}\right]
```
"""

# ╔═╡ d068efed-78d5-4913-a96c-3ec28e5988f7
cm"""
- __Scaled partial pivoting__. Here ``r_k`` is the smallest index such that
```math
\frac{\left|a_{r_k, k}^{(k)}\right|}{s_k}:=\max \left\{\frac{\left|a_{i, k}^{(k)}\right|}{s_k}: k \leq i \leq n\right\}, \quad s_k:=\max _{1 \leq j \leq n}\left|a_{k j}\right|
```
-  __Complete Pivoting__
```math
a_{r_k, s_k}^{(k)}:=\max \left\{\left|a_{i, j}^{(k)}\right|: k \leq i, j \leq n\right\}
```
$(add_space(10))with ``r_k``, ``s_k`` the smallest such indices in case of a tie, is known as. 

- Complete pivoting is known to be more numerically stable than partial pivoting, but requires a lot of search and is seldom used in practice.
"""

# ╔═╡ e0504a1a-d598-4ca7-9023-77d5d57814e0
cm"""
$(ex())Find PLU factorization of
```math
A=\left[\begin{array}{cccc}1 & -1 & 1 & 2 \\ -2 & 1 & 1 & 1 \\ 2 & -1 & 2 & 3 \\ -4 & 1 & 0 & 2\end{array}\right]
```
"""

# ╔═╡ 7090fc3e-97e5-4819-af88-1bb2dd51d7ec
cm"""
$(bbl("Lemma","3.1 (L1U of Leading Principal Submatrices)")) Suppose ``\boldsymbol{A}=\boldsymbol{L} \boldsymbol{U}`` is an L1 ``U`` factorization of ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``. For ``k=1, \ldots, n`` let ``\boldsymbol{A}_{[k]}, \boldsymbol{L}_{[k]}, \boldsymbol{U}_{[k]}`` be the leading principal submatrices of ``\boldsymbol{A}, \boldsymbol{L}, \boldsymbol{U}``, respectively. Then ``\boldsymbol{A}_{[k]}=\boldsymbol{L}_{[k]} \boldsymbol{U}_{[k]}`` is an ``L 1 U`` factorization of ``\boldsymbol{A}_{[k]}`` for ``k=1, \ldots, n``.
"""

# ╔═╡ 802688ba-a8f5-419b-a135-b6525679fb42
cm"""
$(bth("3.4 (LU Theorem)")) A square matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has a unique ``L 1 U`` (LU1, ``L D U)`` factorization if and only if the leading principal submatrices ``\boldsymbol{A}_{[k]}`` of ``\boldsymbol{A}`` are nonsingular for ``k=1, \ldots, n-1``.
"""

# ╔═╡ 81cb8f6c-d82b-41eb-bd06-f4a601954785
cm"""
$(bth("3.5 (Block LU Theorem)")) Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is a block matrix of the form (3.19). Then ``\boldsymbol{A}`` has a unique block ``L U`` factorization (3.20) if and only if the leading principal block submatrices
```math
\boldsymbol{A}_{\{k\}}:=\left[\begin{array}{ccc}
\boldsymbol{A}_{11} & \cdots & \boldsymbol{A}_{1 k} \\
\vdots & & \vdots \\
\boldsymbol{A}_{k 1} & \cdots & \boldsymbol{A}_{k k}
\end{array}\right]
```
are nonsingular for ``k=1, \ldots, m-1``.
"""

# ╔═╡ fb35d4ab-fc81-4618-be09-4a2911ff4566
cm"""
There are special versions of the LU factorization for Hermitian and positive definite. __The most important ones are__
1. the __LDL* factorization__ which is an LDU factorization with ``\boldsymbol{U}=\boldsymbol{L}^*`` and ``\boldsymbol{D}`` a diagonal matrix with real diagonal elements
2. the __LL* factorization__  (called a __Cholesky factorization__) which is an LU factorization with ``\boldsymbol{U}=\boldsymbol{L}^*`` and ``l_{i i}>0`` all ``i``.

A matrix ``\boldsymbol{A}`` having an LDL* factorization must be Hermitian since ``\boldsymbol{D}`` is real so that ``\boldsymbol{A}^*=\left(\boldsymbol{L} \boldsymbol{D} \boldsymbol{L}^*\right)^*=\boldsymbol{L} \boldsymbol{D}^* \boldsymbol{L}^*=\boldsymbol{A}``. The LL* factorization is .

$(ex(1)) (LDL* of ``2 \times 2`` Hermitian Matrix) Let ``a, d \in \mathbb{R}`` and ``b \in \mathbb{C}``. An LDL* factorization of a ``2 \times 2`` Hermitian matrix must satisfy the equations
```math
\left[\begin{array}{ll}
a & \bar{b} \\
b & d
\end{array}\right]=\left[\begin{array}{ll}
1 & 0 \\
l_1 & 1
\end{array}\right]\left[\begin{array}{cc}
d_1 & 0 \\
0 & d_2
\end{array}\right]\left[\begin{array}{ll}
1 & \overline{l_1} \\
0 & 1
\end{array}\right]=\left[\begin{array}{cc}
d_1 & d_1 \overline{l_1} \\
d_1 l_1 & d_1\left|l_1\right|^2+d_2
\end{array}\right]
```
"""

# ╔═╡ 2ad11900-f756-4600-9565-55108bd0e296
cm"""
$(bbl("Lemma","4.1 (LDL* of Leading Principal Sub Matrices)"))
Suppose ``\boldsymbol{A}=\boldsymbol{L} \boldsymbol{D} \boldsymbol{L}^*`` is an ``L D L^*`` factorization of ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``. For ``k=1, \ldots, n`` let ``\boldsymbol{A}_{[k]}, \boldsymbol{L}_{[k]}`` and ``\boldsymbol{D}_{[k]}`` be the leading principal submatrices of ``\boldsymbol{A}, \boldsymbol{L}`` and ``\boldsymbol{D}``, respectively. Then ``\boldsymbol{A}_{[k]}=`` ``\boldsymbol{L}_{[k]} \boldsymbol{D}_{[k]} \boldsymbol{L}_{[k]}^*`` is an ``L D L^*`` factorization of ``\boldsymbol{A}_{[k]}`` for ``k=1, \ldots, n``
"""

# ╔═╡ 3b933d21-8ea9-4790-8588-4662a84481ba
cm"""
$(bth("4.1 (LDL* Theorem)"))
The matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has a unique ``L D L^*`` factorization if and only if ``\boldsymbol{A}=\boldsymbol{A}^*`` and ``\boldsymbol{A}_{[k]}`` is nonsingular for ``k=1, \ldots, n-1``.
"""

# ╔═╡ e90b24e7-1c5a-4190-9af6-34fe939ad47c
cm"""
Given ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``. The function ``f: \mathbb{C}^n \rightarrow \mathbb{R}`` given by
```math
f(\boldsymbol{x})=\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x}=\sum_{i=1}^n \sum_{j=1}^n a_{i j} \bar{x}_i x_j
```
is called a __quadratic form__. 
- Note that ``f`` is real valued if ``\boldsymbol{A}`` is Hermitian. Indeed, ``\overline{f(\boldsymbol{x})}=\overline{\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x}}=\left(\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x}\right)^*=\boldsymbol{x}^* \boldsymbol{A}^* \boldsymbol{x}=f(\boldsymbol{x})``.

$(define("4.1 (Positive Definite Matrix)"))
We say that a matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is

1. __positive definite__ if ``\boldsymbol{A}^*=\boldsymbol{A}`` and ``\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x}>0`` for all nonzero ``\boldsymbol{x} \in \mathbb{C}^n``;
2. __positive semidefinite__ if ``\boldsymbol{A}^*=\boldsymbol{A}`` and ``\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x} \geq 0`` for all ``\boldsymbol{x} \in \mathbb{C}^n``;
3. __negative (semi)definite__ if ``-\boldsymbol{A}`` is positive (semi)definite.
"""


# ╔═╡ 32228ef3-45bb-473c-9a04-21070b475b19
cm"""
$(bbl("Remarks",""))
1. The ``\mathbf{0}`` is positive semidefinite, while the ``\boldsymbol{I}`` is positive definite.
2. The matrix ``\boldsymbol{A}`` is positive definite ``\Longleftrightarrow`` it is positive semidefinite and ``\boldsymbol{x}^* \boldsymbol{A x}=0 \Longrightarrow \boldsymbol{x}=\mathbf{0}``.
3. A positive definite matrix ``\boldsymbol{A}`` is nonsingular. For if ``\boldsymbol{A} \boldsymbol{x}=\mathbf{0}`` then ``\boldsymbol{x}^* \boldsymbol{A} \boldsymbol{x}=0`` and this implies that ``\boldsymbol{x}=\mathbf{0}``.
4. It follows from Lemma 4.6 (below) that a nonsingular positive semidefinite matrix is positive definite.
5. If ``\boldsymbol{A}`` is real then it is enough to show definiteness for real vectors only. Indeed, if ``\boldsymbol{A} \in \mathbb{R}^{n \times n}, \boldsymbol{A}^T=\boldsymbol{A}`` and ``\boldsymbol{x}^T \boldsymbol{A} \boldsymbol{x}>0`` for all nonzero ``\boldsymbol{x} \in \mathbb{R}^n`` then ``\boldsymbol{z}^* \boldsymbol{A} \boldsymbol{z}>0`` for all nonzero ``z \in \mathbb{C}^n``. For if ``\boldsymbol{z}=\boldsymbol{x}+i \boldsymbol{y} \neq \mathbf{0}`` with ``\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n`` then
```math
\begin{aligned}
z^* \boldsymbol{A} z & =(\boldsymbol{x}-i \boldsymbol{y})^T \boldsymbol{A}(\boldsymbol{x}+i \boldsymbol{y})=\boldsymbol{x}^T \boldsymbol{A} \boldsymbol{x}-i \boldsymbol{y}^T \boldsymbol{A} \boldsymbol{x}+i \boldsymbol{x}^T \boldsymbol{A} \boldsymbol{y}-i^2 \boldsymbol{y}^T \boldsymbol{A} \boldsymbol{y} \\
& =\boldsymbol{x}^T \boldsymbol{A} \boldsymbol{x}+\boldsymbol{y}^T \boldsymbol{A} \boldsymbol{y}
\end{aligned}
```
and this is positive since at least one of the real vectors ``\boldsymbol{x}, \boldsymbol{y}`` is nonzero.
"""

# ╔═╡ f7f904ea-2c8b-473f-a37f-53e6ed627a77
cm"""
$(bbl("Lemma","4.2 The Matrix A*A")) The matrix ``\boldsymbol{A}^* \boldsymbol{A}`` is positive semidefinite for any ``m, n \in \mathbb{N}`` and ``\boldsymbol{A} \in \mathbb{C}^{m \times n}``. It is positive definite if and only if ``\boldsymbol{A}`` has linearly independent columns or equivalently rank ``n``.
"""

# ╔═╡ 5225f24c-d76d-4011-8fd3-8c9e16e886c6
cm"""
$(bbl("Lemma","4.3 T is Positive Definite")) The second derivative matrix ``\boldsymbol{T}=`` ``\operatorname{tridiag}(-1,2,-1) \in \mathbb{R}^{n \times n}`` is positive definite.
"""

# ╔═╡ ea4df632-2f6f-4cde-acf7-348b8ebbfb63
cm"""
- Recall that a principal submatrix ``\boldsymbol{B}=\boldsymbol{A}(\boldsymbol{r}, \boldsymbol{r}) \in \mathbb{C}^{k \times k}`` of a matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has elements ``b_{i, j}=a_{r_i, r_j}`` for ``i, j=1, \ldots, k``, where ``1 \leq r_1< \cdots< r_k \leq n``. It is a leading principal submatrix, denoted ``\boldsymbol{A}_{[k]}`` if ``\boldsymbol{r}=[1,2, \ldots, k]^T``. We have
```math
\boldsymbol{A}(\boldsymbol{r}, \boldsymbol{r})=\boldsymbol{X}^* \boldsymbol{A} \boldsymbol{X}, \quad \boldsymbol{X}:=\left[\boldsymbol{e}_{r_1}, \ldots, \boldsymbol{e}_{r_k}\right] \in \mathbb{C}^{n \times k}
```

$(bbl("Lemma","4.4 (Submatrices)"))
Any principal submatrix of a positive (semi)definite matrix is positive (semi)definite.
"""

# ╔═╡ 3a9f7288-bcac-451e-b361-86102a018c3f
cm"""
$(bth("4.2 (LDL* and LL*)")) The following is equivalent for a matrix ``\boldsymbol{A} \in`` ``\mathbb{C}^{n \times n}``.
1. A is positive definite,
2. ``\boldsymbol{A}`` has an ``L D L^*`` factorization with positive diagonal elements in ``\boldsymbol{D}``,
3. A has a Cholesky factorization.

If the Cholesky factorization exists it is unique.
"""

# ╔═╡ dc85a1f4-bff2-466b-98a8-44328d6cf7c1
cm"""
$(ex())
Find the LDL* for the matrix 
```math
A = \begin{bmatrix}
2 &4 &-3\\
4 &14 &-9\\
-3& -9& 12
\end{bmatrix}
```
"""

# ╔═╡ f68a295a-8fe6-47d1-9ca6-9b6199163ed7
cm"""
$(bth("4.3 (Necessary Conditions for Positive (Semi)Definiteness)")) If ``\boldsymbol{A} \in`` ``\mathbb{C}^{n \times n}`` is positive (semi)definite then for all ``i, j`` with ``i \neq j``
1. ``a_{i i}>0,\left(a_{i i} \geq 0\right)``,
2. ``\left|\operatorname{Re}\left(a_{i j}\right)\right|<\left(a_{i i}+a_{j j}\right) / 2,\left(\left|\operatorname{Re}\left(a_{i j}\right)\right| \leq\left(a_{i i}+a_{j j}\right) / 2\right)``,
3. ``\left|a_{i j}\right|<\sqrt{a_{i i} a_{j j}},\left(\left|a_{i j}\right| \leq \sqrt{a_{i i} a_{j j}}\right)``
4. If ``\boldsymbol{A}`` is positive semidefinite and ``a_{i i}=0`` for some ``i`` then ``a_{i j}=a_{j i}=0`` for ``j=1, \ldots, n``.
"""

# ╔═╡ 27e009e2-2e77-4e58-a2d8-50a9efb99e32
cm"""
$(bbl("Lemma","4.5 (Positive Eigenvalues)"))
A matrix is positive (semi)definite if and only if it is Hermitian and all its eigenvalues are positive (nonnegative).
"""

# ╔═╡ 9579da39-7241-4266-918a-0c46b656cbf3
cm"""
$(bbl("Lemma","4.6 (Positive Semidefinite and Nonsingular)")) A matrix is positive definite if and only if it is positive semidefinite and nonsingular.
"""

# ╔═╡ 6e8ecdb1-7d32-429f-9d23-8b79db0cc75d
cm"""
$(bth("4.4 (Positive Definite Characterization)"))
The following statements are equivalent for a matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``.
1. ``\boldsymbol{A}`` is positive definite.
2. ``\boldsymbol{A}`` is Hermitian with only positive eigenvalues.
3. ``\boldsymbol{A}`` is Hermitian and all leading principal submatrices have a positive determinant.
4. ``\boldsymbol{A}=\boldsymbol{B} \boldsymbol{B}^*`` for a nonsingular ``\boldsymbol{B} \in \mathbb{C}^{n \times n}``.
"""

# ╔═╡ a3a98e03-e9ea-424e-972d-f8367cd52642
cm"""
$(define("5.1 (Inner Product)")) An inner product in a complex vector space ``\mathcal{V}`` is a function ``\mathcal{V} \times \mathcal{V} \rightarrow \mathbb{C}`` satisfying for all ``\boldsymbol{x}, \boldsymbol{y}, \boldsymbol{z} \in \mathcal{V}`` and all ``a, b \in \mathbb{C}`` the following conditions: 
1. ``\langle\boldsymbol{x}, \boldsymbol{x}\rangle \geq 0`` with equality if and only if ``\boldsymbol{x}=\mathbf{0}``.
(positivity)
2. ``\langle\boldsymbol{x}, \boldsymbol{y}\rangle=\overline{\langle\boldsymbol{y}, \boldsymbol{x}\rangle}``
(skew symmetry)
3. ``\langle a \boldsymbol{x}+b \boldsymbol{y}, \boldsymbol{z}\rangle=a\langle\boldsymbol{x}, \boldsymbol{z}\rangle+b\langle\boldsymbol{y}, \boldsymbol{z}\rangle``.
(linearity)
"""

# ╔═╡ 35e8b2a4-06d7-4fa4-9b02-63cdffaf1961
cm"""
$(bbl("Remarks",""))
- Recall that  the __standard inner product in ``\mathbb{C}^n``__ is given by
```math
\langle\boldsymbol{x}, \boldsymbol{y}\rangle:=\boldsymbol{y}^* \boldsymbol{x}=\boldsymbol{x}^T \overline{\boldsymbol{y}}=\sum_{j=1}^n x_j \overline{y_j}
```

Note the complex conjugate on ``\boldsymbol{y}``. It is clearly an inner product in ``\mathbb{C}^n``.

- The function
```math
\|\cdot\|: \mathcal{V} \rightarrow \mathbb{R}, \quad x \longmapsto\|x\|:=\sqrt{\langle x, x\rangle}
```
$(add_space(11))is called the inner product norm.

- The inner product norm for the standard inner product is the Euclidian norm ``\|x\|=\|x\|_2=\sqrt{x^* x}``.
"""

# ╔═╡ c0f9964e-b8f8-4eca-8126-64a139e20afc
cm"""
$(bth("5.1 (Cauchy-Schwarz Inequality)")) For any ``\boldsymbol{x}, \boldsymbol{y}`` in a real or complex inner product space
```math
|\langle\boldsymbol{x}, \boldsymbol{y}\rangle| \leq\|\boldsymbol{x}\|\|\boldsymbol{y}\|
```
with equality if and only if ``\boldsymbol{x}`` and ``\boldsymbol{y}`` are linearly dependent.
"""

# ╔═╡ b66c037a-d1cc-42d2-a45c-9f0190b8f28d
cm"""
$(bth("5.2 (Inner Product Norm)")) For all ``\boldsymbol{x}, \boldsymbol{y}`` in an inner product space and all a in ``\mathbb{C}`` we have
1. ``\|\boldsymbol{x}\| \geq 0`` with equality if and only if ``\boldsymbol{x}=\mathbf{0}``.
(positivity)
2. ``\|a \boldsymbol{x}\|=|a|\|x\|``.
(homogeneity)
3. ``\|x+y\| \leq\|x\|+\|y\|``,
(subadditivity)
where ``\|\boldsymbol{x}\|:=\sqrt{\langle\boldsymbol{x}, \boldsymbol{x}\rangle}``.
"""

# ╔═╡ ff626183-36f2-455a-b141-1b9725219f41
cm"""
$(bbl("Remark",""))
In the real case the Cauchy-Schwarz inequality implies that ``-1 \leq \frac{\langle\boldsymbol{x}, \boldsymbol{y}\rangle}{\|\boldsymbol{x}\|\|\boldsymbol{y}\|} \leq 1`` for nonzero ``\boldsymbol{x}`` and ``\boldsymbol{y}``, so there is a unique angle ``\theta`` in ``[0, \pi]`` such that
```math
\cos \theta=\frac{\langle\boldsymbol{x}, \boldsymbol{y}\rangle}{\|\boldsymbol{x}\|\|\boldsymbol{y}\|}
```

This defines the angle between vectors in a real inner product space.
"""

# ╔═╡ e42b6939-e52d-4b9b-adba-02f4549fd08b
cm"""
$(define("5.2 (Orthogonality)"))
Two vectors ``\boldsymbol{x}, \boldsymbol{y}`` in a real or complex inner product space are orthogonal or perpendicular, denoted as ``\boldsymbol{x} \perp \boldsymbol{y}``, if ``\langle\boldsymbol{x}, \boldsymbol{y}\rangle=0``. The vectors are orthonormal if in addition ``\|\boldsymbol{x}\|=\|\boldsymbol{y}\|=1``.
$(ebl())

$(bth("5.3 (Pythagoras)")) For a real or complex inner product space
```math
\|x+y\|^2=\|x\|^2+\|y\|^2, \quad \text { if } \quad x \perp y
```
$(eth())

$(define("5.3 (Orthogonal- and Orthonormal Bases)")) A set of nonzero vectors ``\left\{\boldsymbol{v}_1, \ldots, \boldsymbol{v}_k\right\}`` in a subspace ``\mathcal{S}`` of a real or complex inner product space is an orthogonal basis for ``\mathcal{S}`` if it is a basis for ``\mathcal{S}`` and ``\left\langle\boldsymbol{v}_i, \boldsymbol{v}_j\right\rangle=0`` for ``i \neq j``. It is an orthonormal basis for ``\mathcal{S}`` if it is a basis for ``\mathcal{S}`` and ``\left\langle\boldsymbol{v}_i, \boldsymbol{v}_j\right\rangle=\delta_{i j}`` for all ``i, j``.
"""

# ╔═╡ e53a67cc-c38f-4006-8d20-7bbc17e26f66
cm"""
$(bth("(Gram-Schmidt)")) Let ``\left\{s_1, \ldots, s_k\right\}`` be a basis for a real or complex inner product space ``(\mathcal{S},\langle\cdot, \cdot\rangle)``. Define
```math
\boldsymbol{v}_1:=\boldsymbol{s}_1, \quad \boldsymbol{v}_j:=\boldsymbol{s}_j-\sum_{i=1}^{j-1} \frac{\left\langle\boldsymbol{s}_j, \boldsymbol{v}_i\right\rangle}{\left\langle\boldsymbol{v}_i, \boldsymbol{v}_i\right\rangle} \boldsymbol{v}_i, \quad j=2, \ldots, k
```

Then ``\left\{\boldsymbol{v}_1, \ldots, \boldsymbol{v}_k\right\}`` is an orthogonal basis for ``\mathcal{S}`` and the normalized vectors
```math
\left\{\boldsymbol{u}_1, \ldots, \boldsymbol{u}_k\right\}:=\left\{\frac{\boldsymbol{v}_1}{\left\|\boldsymbol{v}_1\right\|}, \ldots, \frac{\boldsymbol{v}_k}{\left\|\boldsymbol{v}_k\right\|}\right\}
```
form an orthonormal basis for ``\mathcal{S}``.
"""

# ╔═╡ 57a9ebc0-a35f-4900-ba3f-9a7df19ecfbf
cm"""
$(bth("5.5 (Orthogonal Extension of Basis)"))
Suppose ``\mathcal{S} \subset \mathcal{T}`` are finite dimensional subspaces of a vector space ``\mathcal{V}``. An orthogonal basis for ``\mathcal{S}`` can always be extended to an orthogonal basis for ``\mathcal{T}``.
"""

# ╔═╡ 6de81e3d-53b3-458c-a9cf-edcef35e0db3
cm"""
$(bbl("Corollary","5.1 (Extending Orthogonal Vectors to a Basis)")) For ``1 \leq k< n``, a set ``\left\{s_1, \ldots, s_k\right\}`` of nonzero orthogonal vectors in ``\mathbb{R}^n`` or ``\mathbb{C}^n`` can be extended to an orthogonal basis for the whole space.
"""

# ╔═╡ a4a31455-1959-4af0-9b56-ec9dec4c94d5
cm"""
$(bbl("Remarks",""))
- ``\mathcal{S}+\mathcal{T}`` is a vector space.
- Every ``\boldsymbol{v} \in \mathcal{S} \oplus \mathcal{T}`` can be decomposed __uniquely__ in the form ``\boldsymbol{v}=\boldsymbol{s}+\boldsymbol{t}``, where ``\boldsymbol{s} \in \mathcal{S}`` and ``\boldsymbol{t} \in \mathcal{T}``.

- We have
```math
\operatorname{dim}(\mathcal{S} \oplus \mathcal{T})=\operatorname{dim}(\mathcal{S})+\operatorname{dim}(\mathcal{T})
```

- The subspaces ``\mathcal{S}`` and ``\mathcal{T}`` in a direct sum are called __complementary subspaces__.
- An orthogonal sum is a direct sum. For if ``\boldsymbol{v} \in \mathcal{S} \cap \mathcal{T}`` then ``\boldsymbol{v}`` is orthogonal to itself, ``\langle\boldsymbol{v}, \boldsymbol{v}\rangle=0``, which implies that ``\boldsymbol{v}=0``. We often write ``\mathcal{T}:=\mathcal{S}^{\perp}``.
- Suppose ``\boldsymbol{v}=\boldsymbol{s}_0+\boldsymbol{t}_0 \in \mathcal{S} \oplus \mathcal{T}``, where ``\boldsymbol{s}_0 \in \mathcal{S}`` and ``\boldsymbol{t}_0 \in \mathcal{T}``. The vector ``\boldsymbol{s}_0`` is called the __oblique projection__ of ``\boldsymbol{v}`` into ``\mathcal{S}`` along ``\mathcal{T}``. Similarly, the vector ``\boldsymbol{t}_0`` is called the __oblique projection__ of ``v`` into ``\mathcal{T}`` along ``\mathcal{S}``. If ``\mathcal{S} \oplus \mathcal{T}`` is an orthogonal sum then ``\boldsymbol{s}_0`` is called the __orthogonal projection__ of ``\boldsymbol{v}`` into ``\mathcal{S}``. Similarly, ``\boldsymbol{t}_0`` is called the __orthogonal projection__ of ``v`` in ``\mathcal{T}=\mathcal{S}^{\perp}``. The orthogonal projections are illustrated in Fig. 5.2.
"""

# ╔═╡ 489d96e6-14b8-4096-829f-71059ad6d25c
cm"""
$(bth("5.6 (Orthogonal Projection)")) Let ``\mathcal{S}`` and ``\mathcal{T}`` be subspaces of a finite dimensional real or complex vector space ``\mathcal{V}`` with an inner product ``\langle\cdot, \cdot\rangle``. The orthogonal projections ``\boldsymbol{s}_0`` of ``\boldsymbol{v} \in \mathcal{S} \stackrel{\perp}{\oplus} \mathcal{T}`` into ``\mathcal{S}`` and ``\boldsymbol{t}_0`` of ``\boldsymbol{v} \in \mathcal{S} \stackrel{\perp}{\oplus} \mathcal{T}`` into ``\mathcal{T}`` satisfy ``v=s_0+t_0``, and

```math
\left\langle\boldsymbol{s}_0, \boldsymbol{s}\right\rangle=\langle\boldsymbol{v}, \boldsymbol{s}\rangle, \quad \text{for all}\quad \boldsymbol{s} \in \mathcal{S}, \quad\left\langle\boldsymbol{t}_0, \boldsymbol{t}\right\rangle=\langle\boldsymbol{v}, \boldsymbol{t}\rangle, \quad\text{for all}\quad\boldsymbol{t} \in \mathcal{T}.
```
Moreover, if ``\left\{\boldsymbol{v}_1, \ldots, \boldsymbol{v}_k\right\}`` is an orthogonal basis for ``\boldsymbol{S}`` then
```math
\boldsymbol{s}_0=\sum_{i=1}^k \frac{\left\langle\boldsymbol{v}, \boldsymbol{v}_i\right\rangle}{\left\langle\boldsymbol{v}_i, \boldsymbol{v}_i\right\rangle} \boldsymbol{v}_i
```
"""

# ╔═╡ 614c7edf-1133-4e93-a5d0-77747b840ca7
cm"""
$(bbl("Corollary","5.2 (Best Approximation)"))
Let ``\mathcal{S}`` be a subspaces of a finite dimensional real or complex vector space ``\mathcal{V}`` with an inner product ``\langle\cdot, \cdot\rangle`` and corresponding norm ``\|\boldsymbol{v}\|:=\sqrt{\langle\boldsymbol{v}, \boldsymbol{v}\rangle}``. If ``\boldsymbol{s}_0 \in \mathcal{S}`` is the orthogonal projection of ``\boldsymbol{v} \in \mathcal{V}`` then
```math
\left\|v-s_0\right\|<\|v-s\|, \text { for all } s \in \mathcal{S}, \boldsymbol{s} \neq s_0
```
$(ebl())

__Proof__ Let ``s_0 \neq s \in \mathcal{S}`` and ``0 \neq u:=s_0-\boldsymbol{s} \in \mathcal{S}``. It follows from (5.9) that ``\left\langle\boldsymbol{v}-\boldsymbol{s}_0, \boldsymbol{u}\right\rangle=0``. By (5.7) (Pythagoras) we obtain
```math
\|\boldsymbol{v}-\boldsymbol{s}\|^2=\left\|\boldsymbol{v}-\boldsymbol{s}_0+\boldsymbol{u}\right\|^2=\left\|\boldsymbol{v}-\boldsymbol{s}_0\right\|^2+\|\boldsymbol{u}\|^2>\left\|\boldsymbol{v}-\boldsymbol{s}_0\right\|^2
```
"""

# ╔═╡ 310f6fa8-f185-4901-8368-b2b268e40bca
cm"""
$(bbl("Lemma","5.1"))
Let ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` and ``\langle\boldsymbol{x}, \boldsymbol{y}\rangle`` be the standard inner product in ``\mathbb{C}^n``. Then
1. ``\boldsymbol{A}^T=\boldsymbol{A} \Longleftrightarrow\langle\boldsymbol{A} \boldsymbol{x}, \boldsymbol{y}\rangle=\langle\boldsymbol{x}, \overline{\boldsymbol{A}} \boldsymbol{y}\rangle`` for all ``\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n``.
2. ``\boldsymbol{A}^*=\boldsymbol{A} \Longleftrightarrow\langle\boldsymbol{A} \boldsymbol{x}, \boldsymbol{y}\rangle=\langle\boldsymbol{x}, \boldsymbol{A} \boldsymbol{y}\rangle`` for all ``\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n``.
"""

# ╔═╡ ef18b323-485f-48f4-95d0-9063bb6ef2e1
cm"""
- A square matrix ``\boldsymbol{U} \in \mathbb{C}^{n \times n}`` is __unitary__ if 

```math
\boldsymbol{U}^* \boldsymbol{U}=\boldsymbol{I}.
```

- If ``\boldsymbol{U}`` is real then ``\boldsymbol{U}^T \boldsymbol{U}=`` ``\boldsymbol{I}`` and ``\boldsymbol{U}`` is called an orthogonal matrix. Unitary and orthogonal matrices have orthonormal columns.

- If ``\boldsymbol{U}^* \boldsymbol{U}=\boldsymbol{I}`` the matrix ``\boldsymbol{U}`` is nonsingular, ``\boldsymbol{U}^{-1}=\boldsymbol{U}^*`` and therefore ``\boldsymbol{U} \boldsymbol{U}^*=`` ``\boldsymbol{U} \boldsymbol{U}^{-1}=\boldsymbol{I}`` as well. Moreover, both the columns and rows of a unitary matrix of order ``n`` form orthonormal bases for ``\mathbb{C}^n``. We also note that the product of two unitary matrices is unitary. Indeed, if ``\boldsymbol{U}_1^* \boldsymbol{U}_1=\boldsymbol{I}`` and ``\boldsymbol{U}_2^* \boldsymbol{U}_2=\boldsymbol{I}`` then ``\left(\boldsymbol{U}_1 \boldsymbol{U}_2\right)^*\left(\boldsymbol{U}_1 \boldsymbol{U}_2\right)=`` ``\boldsymbol{U}_2^* \boldsymbol{U}_1^* \boldsymbol{U}_1 \boldsymbol{U}_2=\boldsymbol{I}``

$(bth("5.7 (Unitary Matrix)"))
The matrix ``\boldsymbol{U} \in \mathbb{C}^{n \times n}`` is unitary if and only if ``\langle\boldsymbol{U} \boldsymbol{x}, \boldsymbol{U} \boldsymbol{y}\rangle=\langle\boldsymbol{x}, \boldsymbol{y}\rangle`` for all ``\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n``. In particular, if ``\boldsymbol{U}`` is unitary then ``\|\boldsymbol{U} \boldsymbol{x}\|_2=`` ``\|\boldsymbol{x}\|_2`` for all ``\boldsymbol{x} \in \mathbb{C}^n``.
"""

# ╔═╡ 5ae18425-3dde-487d-80f7-127d59a18bbb
cm"""
$(define("5.4 (Householder Transformation)"))
A matrix ``\boldsymbol{H} \in \mathbb{C}^{n \times n}`` of the form
```math
\boldsymbol{H}:=\boldsymbol{I}-\boldsymbol{u} \boldsymbol{u}^*, \text { where } \boldsymbol{u} \in \mathbb{C}^n \text { and } \boldsymbol{u}^* \boldsymbol{u}=2
```
is called a Householder transformation. The name __elementary reflector__ is also used.
"""

# ╔═╡ e28d00ed-710d-4c08-af2e-0f9562d64be2
cm"""
$(bbl("Lemma","5.2"))
Suppose ``\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n`` are two vectors such that ``\|\boldsymbol{x}\|_2=\|\boldsymbol{y}\|_2, \boldsymbol{y}^* \boldsymbol{x}`` is real and ``v:=\boldsymbol{x}-\boldsymbol{y} \neq \mathbf{0}``. Then ``\left(\boldsymbol{I}-2 \frac{v v^*}{v^* v}\right) x=y``.
"""

# ╔═╡ 4eade3b9-589f-492a-9965-03eb74afd493
cm"""
$(post_img("https://www.dropbox.com/scl/fi/cjgji3nv7csvhfmplnxg3/fig5_3.png?rlkey=fbn2wzn62q3ni3gmfxpepq96m&raw=1"))
"""

# ╔═╡ 7627739e-9f51-4197-8ece-ad3e17b0f906
begin
	s52_s = @bind s52s Slider(-4:0.1:5, show_value=true, default=1)
	s52_t = @bind s52t Slider(-4:0.1:5, show_value=true, default=1)
	s52_show =@bind s52show CheckBox(default=false)
	cm"""
	``v_1=`` $(s52_s) 	$(add_space(20))		``v_2=`` $(s52_t)

	Show Zeros of vector ``x`` $(s52_show)
	"""
end

# ╔═╡ a6148859-78c2-4c40-aa0e-0373870a74b8
let
	v=[s52s,s52t]
	u=√2*v/norm(v)
	H(u)=I-2*u*u'/(u'*u)
	x =[-1,2]
	y = H(u)*x
	f(x)=-u[1]*x/u[2]
	p = plot(;frame_style=:origin,xlims=(-5,5),ylims=(-5,5))
	if s52show
	p = plot(p,[0,norm(x)],[0,0],arrow=true,color=:purple,linewidth=2,label="")
	end
	p = plot(p,[0,v[1]],[0,v[2]],arrow=true,color=:red,linewidth=2,label="")
	p = plot(p,[0,x[1]],[0,x[2]],arrow=true,color=:black,linewidth=2,label="")
	p = plot(p,[0,y[1]],[0,y[2]],arrow=true,color=:blue,linewidth=2,label="")
	p = plot(p,z->f(z),c=:green)
	annotate!([(x[1],0.2+x[2],L"x"),(y[1],y[2]-0.2,L"Hx"),(-4.5,f(-4.5),("Mirror",:green)),(v[1],v[2]+0.2,(L"v",:red))])
	if s52show
	annotate!([((norm(x)*1.2,0.2,L"\|x\|e_1"))])
	end
	p
	
end

# ╔═╡ c5077118-527a-448b-bfbc-ad5f00082b7b
cm"""
$(bth("5.8 (Zeros in Vectors)"))
For any ``\boldsymbol{x} \in \mathbb{C}^n`` there is a Householder transformation ``\boldsymbol{H} \in \mathbb{C}^{n \times n}`` such that
```math
\boldsymbol{H} \boldsymbol{x}=a \boldsymbol{e}_1, \quad a=-\rho\|\boldsymbol{x}\|_2, \quad \rho:= \begin{cases}x_1 /\left|x_1\right|, & \text { if } x_1 \neq 0 \\ 1, & \text { otherwise }\end{cases}
```
"""

# ╔═╡ a71dc775-4826-494d-8b83-62274561e6be
cm"""
$(define("Upper Trapezoidal Matrices"))
We say that a matrix ``\boldsymbol{R} \in \mathbb{C}^{m \times n}`` is __upper trapezoidal__, if ``r_{i, j}=0`` for ``j< i`` and ``i=2,3 \ldots, m``. Upper trapezoidal matrices corresponding to ``m< n, m=n``, and ``m>n`` look as follows:
```math
\left[\begin{array}{llll}
x & x & x & x \\
0 & x & x & x \\
0 & 0 & x & x
\end{array}\right], \quad\left[\begin{array}{llll}
x & x & x & x \\
0 & x & x & x \\
0 & 0 & x & x \\
0 & 0 & 0 & x
\end{array}\right], \quad\left[\begin{array}{lll}
x & x & x \\
0 & x & x \\
0 & 0 & x \\
0 & 0 & 0
\end{array}\right] .
```
"""

# ╔═╡ 188ddb61-6cba-4485-83a0-ff37870cebed
cm"""
$(define("QR Decomposition")) Let ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` with ``m, n \in \mathbb{N}``. We say that ``\boldsymbol{A}=\boldsymbol{Q R}`` is a ``\mathbf{Q R}`` decomposition of ``\boldsymbol{A}`` if ``\boldsymbol{Q} \in \mathbb{C}^{m, m}`` is square and unitary and ``\boldsymbol{R} \in \mathbb{C}^{m \times n}`` is upper trapezoidal. If ``m \geq n`` then ``\boldsymbol{R}`` takes the form
```math
\boldsymbol{R}=\left[\begin{array}{c}
\boldsymbol{R}_1 \\
\mathbf{0}_{m-n, n}
\end{array}\right]
```
where ``\boldsymbol{R}_1 \in \mathbb{C}^{n \times n}`` is upper triangular and ``\mathbf{0}_{m-n, n}`` is the zero matrix with ``m-n`` rows and ``n`` columns. For ``m \geq n`` we call ``\boldsymbol{A}=\boldsymbol{Q}_1 \boldsymbol{R}_1`` a ``\mathbf{Q R}`` factorization of ``\boldsymbol{A}`` if ``\boldsymbol{Q}_1 \in \mathbb{C}^{m \times n}`` has orthonormal columns and ``\boldsymbol{R}_1 \in \mathbb{C}^{n \times n}`` is upper triangular.
"""

# ╔═╡ 40820595-a553-4e90-9b78-d6d4b3c473ae
cm"""
$(ex()) Consider
```math
\boldsymbol{A}=\left[\begin{array}{rrr}1 & 3 & 1 \\ 1 & 3 & 7 \\ 1 & -1 & -4 \\ 1 & -1 & 2\end{array}\right]
```
"""

# ╔═╡ 37b9d55e-bc0d-438e-bb99-690b056bd2df
cm"""
$(bth("5.9 (Existence of QR Decomposition)"))
Any matrix ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` with ``m, n \in \mathbb{N}`` has a ``Q R`` decomposition.
"""

# ╔═╡ 74d2a314-2eae-491e-a4ad-dc1184195c00
cm"""
$(bth("5.10 (Uniqueness of QR Factorization)")) If ``m \geq n`` then the ``Q R`` factorization is unique if ``\boldsymbol{A}`` has linearly independent columns and ``\boldsymbol{R}`` has positive diagonal elements.
"""


# ╔═╡ be674a64-cb53-47d0-bdfb-33d72266335d
cm"""
$(ex()) Consider
```math
\boldsymbol{A}=\left[\begin{array}{cc}2 & -1 \\ -1 & 2\end{array}\right]
```
Find ``QR`` descomposition of ``A``.
"""

# ╔═╡ 2e9130ea-3634-45a3-9658-d7b160ed491d
cm"""
> The Gram-Schmidt orthogonalization of the columns of ``\boldsymbol{A}`` can be used to find the QR factorization of ``\boldsymbol{A}``.

$(bth("5.11 (`QR` and Gram-Schmidt)"))
Suppose ``\boldsymbol{A} \in \mathbb{R}^{m \times n}`` has rank ``n`` and let ``\boldsymbol{v}_1, \ldots, \boldsymbol{v}_n`` be the result of applying Gram Schmidt to the columns ``\boldsymbol{a}_1, \ldots, \boldsymbol{a}_n`` of ``\boldsymbol{A}``, i.e.,
```math
\boldsymbol{v}_1=\boldsymbol{a}_1, \quad \boldsymbol{v}_j=\boldsymbol{a}_j-\sum_{i=1}^{j-1} \frac{\boldsymbol{a}_j^T \boldsymbol{v}_i}{\boldsymbol{v}_i^T \boldsymbol{v}_i} \boldsymbol{v}_i, \quad \text { for } j=2, \ldots, n
```

Let
```math
\boldsymbol{Q}_1:=\left[\boldsymbol{q}_1, \ldots, \boldsymbol{q}_n\right], \quad \boldsymbol{q}_j:=\frac{\boldsymbol{v}_j}{\left\|\boldsymbol{v}_i\right\|_2}, \quad j=1, \ldots, n \text { and }
```
```math
\boldsymbol{R}_1:=\left[\begin{array}{cccccc}\left\|\boldsymbol{v}_1\right\|_2 & \boldsymbol{a}_2^T \boldsymbol{q}_1 & \boldsymbol{a}_3^T \boldsymbol{q}_1 & \cdots & \boldsymbol{a}_{n-1}^T \boldsymbol{q}_1 & \boldsymbol{a}_n^T \boldsymbol{q}_1 \\ 0 & \left\|\boldsymbol{v}_2\right\|_2 & \boldsymbol{a}_3^T \boldsymbol{q}_2 & \cdots & \boldsymbol{a}_{n-1}^T \boldsymbol{q}_2 & \boldsymbol{a}_n^T \boldsymbol{q}_2 \\ & 0 & \left\|\boldsymbol{v}_3\right\|_2 & \cdots & \boldsymbol{a}_{n-1}^T \boldsymbol{q}_3 & \boldsymbol{a}_n^T \boldsymbol{q}_3 \\ & & \ddots & \ddots & \vdots & \vdots \\ & & & \ddots & \left\|\boldsymbol{v}_{n-1}\right\|_2 \boldsymbol{a}_n^T \boldsymbol{q}_{n-1} \\ & & & & 0 & \left\|\boldsymbol{v}_n\right\|_2\end{array}\right].
```
Then ``\boldsymbol{A}=\boldsymbol{Q}_1 \boldsymbol{R}_1`` is the unique ``Q R`` factorization of ``\boldsymbol{A}``.
"""

# ╔═╡ f0a58c35-47dc-4666-befc-08c502e6e229
cm"""
$(bbl("Recall",""))
Let ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is a square matrix, ``\lambda \in \mathbb{C}`` and ``\boldsymbol{x} \in \mathbb{C}^n`` then ``(\lambda, \boldsymbol{x})`` is an eigenpair for ``\boldsymbol{A}`` if 
```math 
\boldsymbol{A} \boldsymbol{x}=\lambda \boldsymbol{x}\quad \text{and}\quad  \boldsymbol{x}\quad \text{is nonzero}.
```
- The scalar ``\lambda`` is called an __eigenvalue__ and 
- ``\boldsymbol{x}`` is said to be an __eigenvector__. 
- The set of eigenvalues is called the __spectrum__ of ``\boldsymbol{A}`` and is denoted by ``\sigma(\boldsymbol{A})``. 

For example, ``\sigma(\boldsymbol{I})=\{1, \ldots, 1\}=\{1\}``. 

The eigenvalues are the roots of the __characteristic polynomial__ of ``\boldsymbol{A}`` given for ``\lambda \in \mathbb{C}`` by
```math
\pi_{\boldsymbol{A}}(\lambda)=\operatorname{det}(\boldsymbol{A}-\lambda \boldsymbol{I})
```

The equation ``\operatorname{det}(\boldsymbol{A}-\lambda \boldsymbol{I})=0`` is called the __characteristic equation__ of ``\boldsymbol{A}``. Equivalently the characteristic equation can be written ``\operatorname{det}(\lambda \boldsymbol{I}-\boldsymbol{A})=0``.
"""

# ╔═╡ 3f925be3-06da-4b91-b9c8-1749ba55b3d0
cm"""
$(define(""))
We say that ``\boldsymbol{A}`` is __defective__ if the eigenvectors do not form a basis for ``\mathbb{C}^n`` and __nondefective__ otherwise.
"""

# ╔═╡ 37344680-8fd2-4440-9f03-b7b325c8965c
cm"""
$(bth("6.1 (Distinct Eigenvalues)"))
A matrix with distinct eigenvalues is __nondefective__, i.e., its eigenvectors are linearly independent.
"""

# ╔═╡ ea67c00c-4aed-46ff-a11f-5455553901c0
cm"""
$(ex("Example 6.1","(Defective and Nondefective Matrices)"))
Consider the matrices
```math
\boldsymbol{I}:=\left[\begin{array}{ll}
1 & 0 \\
0 & 1
\end{array}\right], \quad \boldsymbol{J}:=\left[\begin{array}{ll}
1 & 1 \\
0 & 1
\end{array}\right]
```
"""

# ╔═╡ e51e2e47-e389-4afa-b1e0-1ed5bf3768f6
cm"""
$(define("eigenvector expansion"))
If the eigenvectors ``\boldsymbol{x}_1, \ldots, \boldsymbol{x}_n`` form a basis for ``\mathbb{C}^n`` then any ``\boldsymbol{x} \in \mathbb{C}^n`` can be written
```math
\boldsymbol{x}=\sum_{j=1}^n c_j \boldsymbol{x}_j \text { for some scalars } c_1, \ldots, c_n
```
We call this __an eigenvector expansion__ of ``\boldsymbol{x}``. 
"""

# ╔═╡ bfaf9649-e08a-4adf-9794-539277902565
cm"""
$(ex())
```math
A = \left[\begin{array}{cc}2 & -1 \\ -1 & 2\end{array}\right]
```

"""

# ╔═╡ 458a0966-e841-48dc-8e04-6f617fb70b6b
cm"""
$(define("6.1 (Similar Matrices)"))
Two matrices ``\boldsymbol{A}, \boldsymbol{B} \in \mathbb{C}^{n \times n}`` are said to be similar if there is a nonsingular matrix ``\boldsymbol{S} \in \mathbb{C}^{n \times n}`` such that 
```math
\boldsymbol{B}=\boldsymbol{S}^{-1} \boldsymbol{A} \boldsymbol{S}.
```
The transformation ``\boldsymbol{A} \rightarrow \boldsymbol{B}`` is called a similarity transformation. The columns of ``\boldsymbol{S}`` are denoted by ``\boldsymbol{s}_1, \boldsymbol{s}_2, \ldots, \boldsymbol{s}_n``.
"""

# ╔═╡ 823551c5-01bc-4217-82d9-a87bbb7f03fd
cm"""
$(define("Algebraic Multiplicity"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has ``k`` distinct eigenvalues ``\lambda_1, \ldots, \lambda_k`` with multiplicities ``a_1, \ldots, a_k`` so that
```math
\pi_{\boldsymbol{A}}(\lambda):=\operatorname{det}(\boldsymbol{A}-\lambda \boldsymbol{I})=\left(\lambda_1-\lambda\right)^{a_1} \cdots\left(\lambda_k-\lambda\right)^{a_k}, \quad \lambda_i \neq \lambda_j, i \neq j, \sum_{i=1}^k a_i=n
```

The positive integer ``a_i=a\left(\lambda_i\right)=a_A\left(\lambda_i\right)`` is called the __multiplicity__, or more precisely the __algebraic multiplicity__ of the eigenvalue ``\lambda_i``. The multiplicity of an eigenvalue is __simple__ (*double*, *__triple__*) if ``a_i`` is equal to __one__ (*two*, *__three__*).
"""

# ╔═╡ 4f2e0531-99fa-4c48-867b-710c80b42be5
cm"""
$(define("(Geometric Multiplicity)"))
The __geometric multiplicity__ ``g=g(\lambda)=`` ``g_{\boldsymbol{A}}(\lambda)`` of an eigenvalue ``\lambda`` of ``\boldsymbol{A}`` is the dimension of the nullspace ``\mathcal{N}(\boldsymbol{A}-\lambda \boldsymbol{I})`` where 
```math
\mathcal{N}(\boldsymbol{A}-\lambda \boldsymbol{I}):=\left\{\boldsymbol{x} \in \mathbb{C}^n:(\boldsymbol{A}-\lambda \boldsymbol{I}) \boldsymbol{x}=\mathbf{0}\right\}
```
"""

# ╔═╡ 9fdb84a6-2576-4b58-bf2e-6daa75434f1f
cm"""
$(bth("6.2 (Geometric Multiplicity of Similar Matrices)")) Similar matrices have the same eigenvalues with the same algebraic and geometric multiplicities.
"""

# ╔═╡ e4782556-d145-4b88-ba80-359da77b2362
cm"""
$(bth("6.3 (Geometric Multiplicity)"))
We have
1. The geometric multiplicity of an eigenvalue is always bounded above by the algebraic multiplicity of the eigenvalue. ``g_A(\lambda)\leq a_A(\lambda)``.
2. The number of linearly independent eigenvectors of a matrix equals the sum of the geometric multiplicities of the eigenvalues.
3. A matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has ``n`` linearly independent eigenvectors if and only if the algebraic and geometric multiplicity of all eigenvalues are the same.
"""

# ╔═╡ d93f5977-87d5-42cf-b383-8c861a69ccf5
cm"""
$(ex())
```math
A = \left[\begin{array}{lll}6 & 3 & 4 \\ 0 & 6 & 2 \\ 0 & 0 & 7\end{array}\right]
```
"""

# ╔═╡ 2c711e97-7f1f-4c7c-b7c0-3cfb56e912fd
cm"""
$(define("6.3 (Jordan Block)"))
A Jordan block of order ``m``, denoted ``\boldsymbol{J}_m(\boldsymbol{\lambda})`` is an ``m \times m`` matrix of the form
```math
\boldsymbol{J}_m(\lambda):=\left[\begin{array}{cccccc}
\lambda & 1 & 0 & \cdots & 0 & 0 \\
0 & \lambda & 1 & \cdots & 0 & 0 \\
0 & 0 & \lambda & \cdots & 0 & 0 \\
\vdots & & & & \vdots \\
0 & 0 & 0 & \cdots & \lambda & 1 \\
0 & 0 & 0 & \cdots & 0 & \lambda
\end{array}\right]=\lambda \boldsymbol{I}_m+\boldsymbol{E}_m, \quad \boldsymbol{E}_m:=\left[\begin{array}{cccccc}
0 & 1 & 0 & \cdots & 0 & 0 \\
0 & 0 & 1 & \cdots & 0 & 0 \\
0 & 0 & 0 & \cdots & 0 & 0 \\
\vdots & & & & \vdots \\
0 & 0 & 0 & \cdots & 0 & 1 \\
0 & 0 & 0 & \cdots & 0 & 0
\end{array}\right] .
```
"""

# ╔═╡ a657c07f-ab42-4051-8e03-fedb8af62f6b
cm"""
$(bth("6.4 (The Jordan Factorization of a Matrix)"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has ``k`` distinct eigenvalues ``\lambda_1, \ldots, \lambda_k`` of algebraic multiplicities ``a_1, \ldots, a_k`` and geometric multiplicities ``g_1, \ldots, g_k``. There is a nonsingular matrix

``S \in \mathbb{C}^{n \times n}`` such that
```math
\boldsymbol{J}:=\boldsymbol{S}^{-1} \boldsymbol{A} \boldsymbol{S}=\operatorname{diag}\left(\boldsymbol{U}_1, \ldots, \boldsymbol{U}_k\right), \text { with } \boldsymbol{U}_i \in \mathbb{C}^{a_i \times a_i},\tag{*}
```
where each ``\boldsymbol{U}_i`` is block diagonal having ``g_i`` Jordan blocks along the diagonal
```math
\boldsymbol{U}_i=\operatorname{diag}\left(\boldsymbol{J}_{m_{i, 1}}\left(\lambda_i\right), \ldots, \boldsymbol{J}_{m_{i, g_i}}\left(\lambda_i\right)\right)
```

Here ``m_{i, 1}, \ldots, m_{i, g_i}`` are positive integers and they are unique if they are ordered so that ``m_{i, 1} \geq m_{i, 2} \geq \cdots \geq m_{i, g_i}``. Moreover, ``a_i=\sum_{j=1}^{g_i} m_{i, j}`` for all ``i``.
$(eth())
$(bbl("Remarks"))
We note that
1. The matrices ``\boldsymbol{S}`` and ``\boldsymbol{J}`` in (*) are called Jordan factors. We also call ``\boldsymbol{J}`` the Jordan factorization of ``\boldsymbol{A}``.
2. The columns of ``S`` are called principal vectors or generalized eigenvectors. They satisfy the matrix equation ``\boldsymbol{A} \boldsymbol{S}=\boldsymbol{S J}``.
3. Each ``\boldsymbol{U}_i`` is upper triangular with the eigenvalue ``\lambda_i`` on the diagonal and consists of ``g_i`` Jordan blocks. These Jordan blocks can be taken in any order and it is customary to refer to any such block diagonal matrix as the Jordan factorization of ``\boldsymbol{A}``.
"""

# ╔═╡ b5c411da-2629-4942-944b-5f3f4245b74e
cm"""
$(ex())
```math
A = \left[\begin{array}{lll}6 & 3 & 4 \\ 0 & 6 & 2 \\ 0 & 0 & 7\end{array}\right]
```
Find Jordan factors ``J`` and ``S`` of ``A``.
"""

# ╔═╡ db3892ca-816a-4982-adb0-ef09109a99a0
cm"""
$(ex())
```math
A=\left[\begin{array}{rrrrr}2 & 1 & -1 & 8 & -3 \\ 0 & 2 & 0 & 7 & 5 \\ 0 & 0 & 2 & 7 & 5 \\ 0 & 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 0 & 2\end{array}\right]
```
"""

# ╔═╡ ce0f671e-395d-4cd4-a95b-5ec67e7c1fd0
cm"""
$(ex())
Find Jodan factors of A
```math
A=\left[\begin{array}{llll}2 & 0 & 0 & 0 \\ 0 & 2 & 0 & 0 \\ 0 & 0 & 2 & 1 \\ 1 & 0 & 0 & 2\end{array}\right]
```
"""

# ╔═╡ 4830028b-2e80-48b1-8dc0-79597852b828
cm"""
$(bth("6.5 (Schur Factorization)"))
For each ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` there exists a unitary matrix ``\boldsymbol{U} \in \mathbb{C}^{n \times n}`` such that ``\boldsymbol{R}:=\boldsymbol{U}^* \boldsymbol{A} \boldsymbol{U}`` is upper triangular.

The matrices ``\boldsymbol{U}`` and ``\boldsymbol{R}`` in the __Schur factorization__ are called __Schur factors__. 

We call ``\boldsymbol{A}=\boldsymbol{U} \boldsymbol{R} \boldsymbol{U}^*`` the Schur factorization of ``\boldsymbol{A}``.
"""

# ╔═╡ 70686047-be6a-4d77-8cc3-dcb38b7bc378
cm"""
$(define("Normal Matrix"))
A matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is normal if ``\boldsymbol{A}^* \boldsymbol{A}=\boldsymbol{A} \boldsymbol{A}^*``. In this section we show that a matrix has orthogonal eigenvectors if and only if it is normal.

Examples of normal matrices are
1. ``A^*=A``,
(Hermitian)
2. ``\boldsymbol{A}^*=-\boldsymbol{A}``,
(Skew-Hermitian)
3. ``\boldsymbol{A}^*=\boldsymbol{A}^{-1}``,
(Unitary)
4. ``\boldsymbol{A}=\operatorname{diag}\left(d_1, \ldots, d_n\right)``.
(Diagonal)
"""

# ╔═╡ afa0532c-de60-44d4-a166-befcdc1959e4
cm"""
$(bbl("","(Spectral Theorem for Normal Matrices)"))
A matrix ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is normal if and only if there exists a unitary matrix ``\boldsymbol{U} \in \mathbb{C}^{n \times n}`` (that is ``\boldsymbol{U}^*\boldsymbol{U} =\boldsymbol{I} ``) such that 
```math
\boldsymbol{U}^* \boldsymbol{A} \boldsymbol{U}=\boldsymbol{D} \quad \text{is diagonal.}
```
If ``\boldsymbol{D}=\operatorname{diag}\left(\lambda_1, \ldots, \lambda_n\right)`` and ``\boldsymbol{U}=\left[\boldsymbol{u}_1, \ldots, \boldsymbol{u}_n\right]`` then ``\left(\lambda_j, \boldsymbol{u}_j\right), j=`` ``1, \ldots, n`` are orthonormal eigenpairs for ``\boldsymbol{A}``.
"""

# ╔═╡ b23fc699-d1e8-4c83-8e61-06dc4ab8ea6a
cm"""
$(define("The singular value decomposition (SVD)"))
The singular value decomposition (SVD) is a __decomposition__ of a matrix in the form 
```math
\boldsymbol{A}=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*,
``` 
where 
- ``\boldsymbol{U}`` and ``\boldsymbol{V}`` are unitary, and 
- ``\boldsymbol{\Sigma}`` is a nonnegative diagonal matrix, i.e., ``\Sigma_{i j}=0`` for all ``i \neq j`` and ``\Sigma_{i i} \geq 0`` for all ``i``. 

The diagonal elements ``\sigma_i:=\Sigma_{i i}`` are __called singular values__, while the columns of ``\boldsymbol{U}`` and ``\boldsymbol{V}`` are called __singular vectors__. 

To be a singular value decomposition the singular values should be ordered, i.e., 
```math
\sigma_i \geq \sigma_{i+1} \quad \text{for all}\quad  i.
```
"""

# ╔═╡ a2dc9322-049b-4437-a29e-1cd23f816059
cm"""
$(bbl("Remark",""))
The singular values of an ``m \times n`` matrix ``A`` are the positive square roots of the nonzero eigenvalues of the ``n \times n`` symmetric matrix ``A^t A``.
"""

# ╔═╡ 8e1f15f2-8273-400d-97a1-86cd7d138a64
cm"""
$(bth("7.1 "))
Suppose ``m, n \in \mathbb{N}`` and ``\boldsymbol{A} \in \mathbb{C}^{m \times n}``.
1. The matrices ``\boldsymbol{A}^* \boldsymbol{A} \in \mathbb{C}^{n \times n}`` and ``\boldsymbol{A} \boldsymbol{A}^* \in \mathbb{C}^{m \times m}`` have the same nonzero eigenvalues with the same algebraic multiplicities. Moreover the extra eigenvalues of the larger matrix are all zero.
2. The matrices ``\boldsymbol{A}^* \boldsymbol{A}`` and ``\boldsymbol{A} \boldsymbol{A}^*`` are Hermitian with nonnegative eigenvalues.
3. Let ``\left(\lambda_j, \boldsymbol{v}_j\right)`` be orthonormal eigenpairs for ``\boldsymbol{A}^* \boldsymbol{A}`` with
```math
\lambda_1 \geq \cdots \geq \lambda_r>0=\lambda_{r+1}=\cdots=\lambda_n .
```
$(add_space(10))Then ``\left\{\boldsymbol{A} \boldsymbol{v}_1, \ldots, \boldsymbol{A} \boldsymbol{v}_r\right\}`` is an orthogonal basis for the column space 
```math 
\mathcal{R}(\boldsymbol{A}):=\left\{\boldsymbol{A} \boldsymbol{y} \in \mathbb{C}^m: \boldsymbol{y} \in \mathbb{C}^n\right\}
``` 
$(add_space(10))and ``\left\{\boldsymbol{v}_{r+1}, \ldots, \boldsymbol{v}_n\right\}`` is an orthonormal basis for the nullspace 
```math
\mathcal{N}(\boldsymbol{A}):=\left\{\boldsymbol{y} \in \mathbb{C}^n: \boldsymbol{A} \boldsymbol{y}=\mathbf{0}\right\}.
```
4. Let ``\left(\lambda_j, \boldsymbol{u}_j\right)`` be orthonormal eigenpairs for ``\boldsymbol{A} \boldsymbol{A}^*``. If ``\lambda_j>0, j=1, \ldots, r`` and ``\lambda_j=0, j=r+1, \ldots, m`` then ``\left\{\boldsymbol{A}^* \boldsymbol{u}_1, \ldots, \boldsymbol{A}^* \boldsymbol{u}_r\right\}`` is an orthogonal basis for the column space ``\mathcal{R}\left(\boldsymbol{A}^*\right)`` and ``\left\{\boldsymbol{u}_{r+1}, \ldots, \boldsymbol{u}_m\right\}`` is an orthonormal basis for the nullspace ``\mathcal{N}\left(\boldsymbol{A}^*\right)``.
5. The rank of ``\boldsymbol{A}`` equals the number of positive eigenvalues of ``\boldsymbol{A}^* \boldsymbol{A}`` and ``\boldsymbol{A} \boldsymbol{A}^*``.
"""

# ╔═╡ e1d0a740-2c6b-4b72-ab3e-0c0fde78ed68
cm"""
$(bth("7.2 (Existence of SVD)")) Suppose for ``m, n, r \in \mathbb{N}`` that ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` has rank ``r``, and that ``\left(\lambda_j, \boldsymbol{v}_j\right)`` are orthonormal eigenpairs for ``\boldsymbol{A}^* \boldsymbol{A}`` with ``\lambda_1 \geq \cdots \geq`` ``\lambda_r>0=\lambda_{r+1}=\cdots=\lambda_n``. Define
1. ``\boldsymbol{V}:=\left[\boldsymbol{v}_1, \ldots, \boldsymbol{v}_n\right] \in \mathbb{C}^{n \times n}``,
2. ``\Sigma \in \mathbb{R}^{m \times n}`` is a diagonal matrix with diagonal elements ``\sigma_j:=\sqrt{\lambda_j}`` for ``j=`` ``1, \ldots, \min (m, n)``,
3. ``\boldsymbol{U}:=\left[\boldsymbol{u}_1, \ldots, \boldsymbol{u}_m\right] \in \mathbb{C}^{m \times m}``, where ``\boldsymbol{u}_j=\sigma_j^{-1} \boldsymbol{A} \boldsymbol{v}_j`` for ``j=1, \ldots, r`` and ``\boldsymbol{u}_{r+1}, \ldots, \boldsymbol{u}_m`` is any extension of ``\boldsymbol{u}_1, \ldots, \boldsymbol{u}_r`` to an orthonormal basis ``\boldsymbol{u}_1, \ldots, \boldsymbol{u}_m`` for ``\mathbb{C}^m``.
Then ``\boldsymbol{A}=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*`` is a singular value decomposition of ``\boldsymbol{A}``.
"""

# ╔═╡ b1b78dd1-112f-42b9-9d40-66123d64feac
cm"""
$(ex())
Determine the singular values of the ``2 \times 3`` matrix
```math
A=\left[\begin{array}{lll}
1 & 1 & 0 \\
0 & 0 & 1 \\
\end{array}\right]
```
"""

# ╔═╡ 3013da2b-e879-4730-a46f-ab3a94e6cbb8
cm"""
$(ex())
Determine the singular values of the ``5 \times 3`` matrix
```math
A=\left[\begin{array}{lll}
1 & 0 & 1 \\
0 & 1 & 0 \\
0 & 1 & 1 \\
0 & 1 & 0 \\
1 & 1 & 0
\end{array}\right]
```
"""

# ╔═╡ db97a178-677b-48f7-b15d-fbe64d374128
cm"""
$(bth("7.3 (Singular Vectors and Orthonormal Bases)"))
For positive integers ``m, n`` let ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` have rank ``r`` and a singular value decomposition ``\boldsymbol{A}=`` ``\left[\boldsymbol{u}_1, \ldots, \boldsymbol{u}_m\right] \boldsymbol{\Sigma}\left[\boldsymbol{v}_1, \ldots, \boldsymbol{v}_n\right]^*=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^*``. Then the singular vectors satisfy
```math
\begin{aligned}
\boldsymbol{A} \boldsymbol{v}_i & =\sigma_i \boldsymbol{u}_i, i=1, \ldots, r, & & \boldsymbol{A} \boldsymbol{v}_i=0, i=r+1, \ldots, n \\
\boldsymbol{A}^* \boldsymbol{u}_i & =\sigma_i \boldsymbol{v}_i, i=1, \ldots, r, & & \boldsymbol{A}^* \boldsymbol{u}_i=0, i=r+1, \ldots, m
\end{aligned}
```

Moreover,
1. ``\left\{\boldsymbol{u}_1, \ldots, \boldsymbol{u}_r\right\}`` is an orthonormal basis for ``\mathcal{R}(\boldsymbol{A})``,
2. ``\left\{\boldsymbol{u}_{r+1}, \ldots, \boldsymbol{u}_m\right\}`` is an orthonormal basis for ``\mathcal{N}\left(\boldsymbol{A}^*\right)``,
3. ``\left\{\boldsymbol{v}_1, \ldots, \boldsymbol{v}_r\right\}`` is an orthonormal basis for ``\mathcal{R}\left(\boldsymbol{A}^*\right)``,
4. ``\left\{\boldsymbol{v}_{r+1}, \ldots, \boldsymbol{v}_n\right\}`` is an orthonormal basis for ``\mathcal{N}(\boldsymbol{A})``.
"""

# ╔═╡ 6cbcc482-d188-4732-8926-2cae312e2788
cm"""
$(ex())
```math
\boldsymbol{A}:=\frac{1}{25}\left[\begin{array}{ll}11 & 48 \\ 48 & 39\end{array}\right]
```
"""

# ╔═╡ 56f09d0a-67c5-4109-bde0-b66296986a1c
cm"""
$(define("8.1 (Vector Norm)")) 
A __(vector) norm__ in a real (resp. complex) vector space ``\mathcal{V}`` is a function ``\|\cdot\|: \mathcal{V} \rightarrow \mathbb{R}`` that satisfies for all ``\boldsymbol{x}, \boldsymbol{y}`` in ``\mathcal{V}`` and all ``a`` in ``\mathbb{R}`` (resp. ``\mathbb{C}`` )
1. ``\|\boldsymbol{x}\| \geq 0`` with equality if and only if ``\boldsymbol{x}=\mathbf{0}``.
(positivity)
2. ``\|a x\|=|a|\|x\|``.
(homogeneity)
3. ``\|\boldsymbol{x}+\boldsymbol{y}\| \leq\|\boldsymbol{x}\|+\|\boldsymbol{y}\|``.
(subadditivity)

The triple ``(\mathcal{V}, \mathbb{R},\|\cdot\|)(`` resp. ``(\mathcal{V}, \mathbb{C},\|\cdot\|))`` is called a __normed vector space__ and the inequality 3 . is called the triangle inequality.
"""

# ╔═╡ 4f107e76-56a3-4b87-a76b-ec2179b082df

cm"""
$(define("8.2 (Vector p-Norms)"))
We define for ``p \geq 1`` and ``\boldsymbol{x} \in \mathbb{R}^n`` or ``\boldsymbol{x} \in \mathbb{C}^n`` the ``p``-norms by
```math
\begin{aligned}
\|\boldsymbol{x}\|_p & :=\left(\sum_{j=1}^n\left|x_j\right|^p\right)^{1 / p} \\
\|\boldsymbol{x}\|_{\infty} & :=\max _{1 \leq j \leq n}\left|x_j\right|
\end{aligned}
```
"""

# ╔═╡ 47337c99-d8e1-4488-9959-77963d555e7d
cm"""
$(bbl("Remarks",""))
1. It can be shown that the ``p``-norms are vector norms for ``1 \leq p \leq \infty``.

2. The triangle inequality called Minkowski's inequality
```math
\|\boldsymbol{x}+\boldsymbol{y}\|_p \leq\|\boldsymbol{x}\|_p+\|\boldsymbol{y}\|_p.
```
3. To prove it one first establishes Hölder's inequality
```math
\sum_{j=1}^n\left|x_j y_j\right| \leq\|\boldsymbol{x}\|_p\|\boldsymbol{y}\|_q, \quad \frac{1}{p}+\frac{1}{q}=1, \quad \boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n
```
$(add_space(10))The relation ``\frac{1}{p}+\frac{1}{q}=1`` means that if ``p=1`` then ``q=\infty`` and vice versa. The Hölder's inequality is the same as the Cauchy-Schwarz inequality for the Euclidian norm ``p=2``.

4. The infinity norm is related to the other ``p``-norms by
```math
\lim _{p \rightarrow \infty}\|\boldsymbol{x}\|_p=\|\boldsymbol{x}\|_{\infty} \text { for all } \boldsymbol{x} \in \mathbb{C}^n\tag{*}
```
5. The equation (*) clearly holds for ``\boldsymbol{x}=\mathbf{0}``. For ``\boldsymbol{x} \neq \mathbf{0}`` we write
```math
\|\boldsymbol{x}\|_p:=\|\boldsymbol{x}\|_{\infty}\left(\sum_{j=1}^n\left(\frac{\left|x_j\right|}{\|\boldsymbol{x}\|_{\infty}}\right)^p\right)^{1 / p} .
```
$(add_space(10))Now each term in the sum is not greater than one and at least one term is equal to one, and we obtain
```math
\|\boldsymbol{x}\|_{\infty} \leq\|\boldsymbol{x}\|_p \leq n^{1 / p}\|\boldsymbol{x}\|_{\infty}, \quad p \geq 1\tag{**}
```

$(add_space(10))Since ``\lim _{p \rightarrow \infty} n^{1 / p}=1`` for any fixed ``n \in \mathbb{N}``, we see that (*) follows.

6. We can show the following generalization of inequality (**)
```math
\|\boldsymbol{x}\|_{p^{\prime}} \leq\|\boldsymbol{x}\|_p \leq n^{1 / p-1 / p^{\prime}}\|\boldsymbol{x}\|_{p^{\prime}}, \quad \boldsymbol{x} \in \mathbb{C}^n, \quad 1 \leq p \leq p^{\prime} \leq \infty .
```

"""

# ╔═╡ 53af53ab-a701-47cc-8979-14d6bef5e074
cm"""
$(define("8.3 (Equivalent Norms)"))
We say that two norms ``\|\cdot\|`` and ``\|\cdot\|^{\prime}`` on ``\mathcal{V}`` are equivalent if there are positive constants ``m`` and ``M`` such that for all vectors ``\boldsymbol{x} \in \mathcal{V}`` we have
```math
m\|\boldsymbol{x}\|^{\prime} \leq\|\boldsymbol{x}\| \leq M\|\boldsymbol{x}\|^{\prime} .
```
$(ebl())
$(bbl("Remark",""))
By (**) the ``p`` - and ``\infty``-norms are equivalent for any ``p \geq 1``. This result is generalized in the following theorem.
$(ebl())
$(bth(" 8.1 (Basic Properties of Vector Norms)"))
The following holds for a normed vector space ``(\mathcal{V}, \mathbb{C},\|\cdot\|)``.
1. ``\|\boldsymbol{x}-\boldsymbol{y}\| \geq|\|\boldsymbol{x}\|-\|\boldsymbol{y}\||``, for all ``\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n \quad`` (inverse triangle inequality).
2. The vector norm is a continuous function ``\mathcal{V} \rightarrow \mathbb{R}``.
3. All vector norms on ``\mathcal{V}`` are equivalent provided ``\mathcal{V}`` is finite dimensional.
"""

# ╔═╡ d72010fc-cd62-4313-a3b2-aea4f5b277c2
cm"""
$(bth("8.2 (Matrix Norm Equivalence)"))
All matrix norms on ``\mathbb{C}^{m \times n}`` are equivalent. Thus, if ``\|\cdot\|`` and ``\|\cdot\|^{\prime}`` are two matrix norms on ``\mathbb{C}^{m \times n}`` then there are positive constants ``\mu`` and ``M`` such that
```math
\mu\|\boldsymbol{A}\| \leq\|\boldsymbol{A}\|^{\prime} \leq M\|\boldsymbol{A}\|
```
holds for all ``\boldsymbol{A} \in \mathbb{C}^{m \times n}``. Moreover, a matrix norm is a continuous function.
"""

# ╔═╡ 4da69a7c-1ece-485d-a353-cb944c900b26
cm"""
$(bbl("Lemma","(Frobenius Norm Properties)")) 
For any ``m, n \in \mathbb{N}`` and any matrix ``\boldsymbol{A} \in`` ``\mathbb{C}^{m \times n}``
1. ``\left\|\boldsymbol{A}^*\right\|_F=\|\boldsymbol{A}\|_F``,
2. ``\|\boldsymbol{A}\|_F^2=\sum_{j=1}^n\left\|\boldsymbol{a}_{: j}\right\|_2^2``,
3. ``\|\boldsymbol{U} \boldsymbol{A}\|_F=\|\boldsymbol{A} \boldsymbol{V}\|_F=\|\boldsymbol{A}\|_F`` for any unitary matrices ``\boldsymbol{U} \in \mathbb{C}^{m \times m}`` and ``\boldsymbol{V} \in`` ``\mathbb{C}^{n \times n}``,
4. ``\|\boldsymbol{A} \boldsymbol{B}\|_F \leq\|\boldsymbol{A}\|_F\|\boldsymbol{B}\|_F`` for any ``\boldsymbol{B} \in \mathbb{C}^{n, k}, \quad k \in \mathbb{N}``,
5. ``\|\boldsymbol{A} \boldsymbol{x}\|_2 \leq\|\boldsymbol{A}\|_F\|\boldsymbol{x}\|_2``, for all ``\boldsymbol{x} \in \mathbb{C}^n``.
"""

# ╔═╡ 9a9975ef-2ec1-41af-bd6f-0a76756c100b
cm"""
$(bth("(Frobenius Norm and Singular Values)"))
We have ``\|\boldsymbol{A}\|_F=`` ``\sqrt{\sigma_1^2+\cdots+\sigma_n^2}``, where ``\sigma_1, \ldots, \sigma_n`` are the singular values of ``\boldsymbol{A}``.

__Proof__ Using previous Lemma we find
```math
\|\boldsymbol{A}\|_F \stackrel{\text { 3. }}{=}\left\|\boldsymbol{U}^* \boldsymbol{A} \boldsymbol{V}\right\|_F=\|\boldsymbol{\Sigma}\|_F=\sqrt{\sigma_1^2+\cdots+\sigma_n^2} .
```
"""

# ╔═╡ 263b843b-e75c-45a4-98a9-4f9111858ae6
cm"""
$(define("8.4 (Consistent Matrix Norms)"))
A matrix norm is called consistent on ``\mathbb{C}^{n \times n}`` if

```math
\|\boldsymbol{A} \boldsymbol{B}\| \leq\|\boldsymbol{A}\|\|\boldsymbol{B}\|\quad
\text{(submultiplicativity)}\quad
\text{holds for all } \boldsymbol{A}, \boldsymbol{B} \in \mathbb{C}^{n \times n}.
```
"""

# ╔═╡ fba0bff1-e7a9-4793-94bd-968cd82fd926
cm"""
$(bbl("Remark",""))
For a consistent matrix norm on ``\mathbb{C}^{n \times n}`` we have the inequality 
```math
\left\|\boldsymbol{A}^k\right\| \leq\|\boldsymbol{A}\|^k \text{ for } \boldsymbol{A} \in \mathbb{C}^{n \times n} \text{ and } k \in \mathbb{N}.\tag{🌟}
```
"""

# ╔═╡ 048fc935-3507-496b-b00d-09078784fe06
cm"""
$(define("8.5 (Subordinate Matrix Norms)"))
Suppose ``m, n \in \mathbb{N}`` are given, let ``\| \|`` on ``\mathbb{C}^m`` and ``\left\|\|_\beta\right.`` on ``\mathbb{C}^n`` be vector norms, and let ``\| \|`` be a matrix norm on ``\mathbb{C}^{m \times n}``. 
We say that the matrix norm ``\|\|`` is __subordinate__ to the vector norms ``\| \|`` and ``\left\|\|_\beta\right.`` if 
```math
\|\boldsymbol{A x}\| \leq\|\boldsymbol{A}\|\|\boldsymbol{x}\|_\beta\quad \text{ for all}\quad \boldsymbol{A} \in \mathbb{C}^{m \times n}\quad \text{and all } \boldsymbol{x} \in \mathbb{C}^n.
```
"""

# ╔═╡ 1a8ed4b8-74ab-4775-b01b-24b637d9288c
cm"""
$(bbl("Remark",""))

The Frobenius  norm is __subordinate__ to the Euclidian vector norm

``\color{white}{.}``
"""

# ╔═╡ 40940718-bf05-4f6e-8244-7716c77ca10c
cm"""
$(bbl("Proposition"," 8.1"))
For ``m, n \in \mathbb{N}, \boldsymbol{A} \in \mathbb{C}^{m \times n}``, all ``\boldsymbol{x} \in \mathbb{C}^n`` and any consistent matrix norm || ||
```math
\|\boldsymbol{A x}\| \leq\|\boldsymbol{A}\|\|\boldsymbol{x}\|,
```
i.e., a consistent matrix norm is subordinate to itself. Moreover, the matrix power bound (🌟) holds for all square matrices ``\boldsymbol{A} \in \mathbb{C}^{n \times n}``.

"""

# ╔═╡ 46991e0d-8a9c-487c-a157-72690cf2b489
cm"""
$(define("8.6 (Operator Norm)"))
Let ``\left\|\|\right.`` be a vector norm defined on ``\mathbb{C}^n`` for all ``n \in \mathbb{N}``. For given ``m, n \in \mathbb{N}`` and ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` we define
```math
\|\boldsymbol{A}\|:=\max _{\boldsymbol{x} \neq 0} \frac{\|\boldsymbol{A} \boldsymbol{x}\|}{\|\boldsymbol{x}\|} .\tag{🌻}
```

We call this the __operator norm__ corresponding to the vector norm ``\| \|``.
"""

# ╔═╡ 89819051-5d7b-47f7-ad1c-94362c4679fe
cm"""
$(bbl("Lemma", "8.1 (The Operator Norm Is a Consistent Matrix Norm)")) If ``|| ||`` is vector norm defined on ``\mathbb{C}^n`` for all ``n \in \mathbb{N}``, then the operator norm given by (🌻) is a consistent matrix norm. Moreover, ``\|\boldsymbol{I}\|=1``.
"""

# ╔═╡ a82fb7fd-861d-4843-8cc6-18a33d9a41b4
cm"""
We define for any ``1 \leq p \leq \infty``
```math
\|\boldsymbol{A}\|_p:=\max _{\boldsymbol{x} \neq 0} \frac{\|\boldsymbol{A} \boldsymbol{x}\|_p}{\|\boldsymbol{x}\|_p}=\max _{\|\boldsymbol{y}\|_p=1}\|\boldsymbol{A} \boldsymbol{y}\|_p .
```

$(bth("8.3 (One-Two-Inf-Norms)"))
For ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` we have
```math
\begin{array}{ll}
\|\boldsymbol{A}\|_1:=\max _{1 \leq j \leq n}\left\|\boldsymbol{A} \boldsymbol{e}_j\right\|_1=\max _{1 \leq j \leq n} \sum_{k=1}^m\left|a_{k, j}\right|, \quad \text { (max column sum) } \\
\|\boldsymbol{A}\|_2:=\sigma_1, \quad \quad \text { (largest singular value of } \boldsymbol{A} \text { ) } \\
\|\boldsymbol{A}\|_{\infty}=\max _{1 \leq k \leq m}\left\|\boldsymbol{e}_k^T \boldsymbol{A}\right\|_1=\max _{1 \leq k \leq m} \sum_{j=1}^n\left|a_{k, j}\right| . \quad \quad \text { (max row sum) }
\end{array}
```
 The __two-norm__ ``\|A\|_2`` is also called the __spectral norm__ of ``A``
"""

# ╔═╡ 01efad8f-6aea-45ee-b084-9b68b35fb61b
cm"""
$(ex()) Let 
```math
\boldsymbol{A}:=\frac{1}{15}\left[\begin{array}{ccc}14 & 4 & 16 \\ 2 & 22 & 13\end{array}\right].
```
Then find (i) ``\|A\|_1``, (ii) ``\|A\|_2``, (iii) ``\|A\|_{\infty}`` and (iv) ``\|A\|_F``.
"""

# ╔═╡ 2b4ca2be-25dd-4f9c-a2c2-c4845c0f048a
cm"""
$(bth("8.4 (Spectral Norm)"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` has singular values ``\sigma_1 \geq \sigma_2 \geq`` ``\cdots \geq \sigma_n`` and eigenvalues ``\left|\lambda_1\right| \geq\left|\lambda_2\right| \geq \cdots \geq\left|\lambda_n\right|``. Then
```math 
\begin{array}{lcl}
\|\boldsymbol{A}\|_2&=&\sigma_1 \quad\text{and}\quad \left\|\boldsymbol{A}^{-1}\right\|_2=\frac{1}{\sigma_n},\\
\\
\|\boldsymbol{A}\|_2&=&\lambda_1 \quad\text{and}\quad \left\|\boldsymbol{A}^{-1}\right\|_2=\frac{1}{\lambda_n},\quad  \text{if } \boldsymbol{A} \text{  is positive definite,} \\
\\
\|\boldsymbol{A}\|_2&=&\left|\lambda_1\right| \quad\text{and}\quad \left\|\boldsymbol{A}^{-1}\right\|_2=\frac{1}{\left|\lambda_n\right|},\quad \text{if } \boldsymbol{A} \text{  is normal}.
\end{array}
```
For the norms of ``\boldsymbol{A}^{-1}`` we assume that ``\boldsymbol{A}`` is nonsingular.
"""

# ╔═╡ 56f568bb-f68e-4c92-b1fa-89ed1a422c8a
cm"""
$(bth("8.5 (Spectral Norm Bound)"))
For any ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` we have ``\|\boldsymbol{A}\|_2^2 \leq`` ``\|\boldsymbol{A}\|_1\|\boldsymbol{A}\|_{\infty}``.
"""

# ╔═╡ bbf4261d-7747-47e4-b69d-1d293c6172c0
cm"""
$(define("(Unitary Invariant Norm)"))
A matrix norm ``\left\|\|\right.`` on ``\mathbb{C}^{m \times n}`` is called unitary invariant if ``\|\boldsymbol{U} \boldsymbol{A} \boldsymbol{V}\|=\|\boldsymbol{A}\|`` for any ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` and any unitary matrices ``\boldsymbol{U} \in \mathbb{C}^{m \times m}`` and ``\boldsymbol{V} \in \mathbb{C}^{n \times n}``.
"""

# ╔═╡ 076b87bb-261d-4dca-a520-67f760a22d84
cm"""
$(bth("8.6 (Unitary Invariant Norms)"))
The Frobenius norm and the spectral norm are unitary invariant. Moreover,
```math
\left\|\boldsymbol{A}^*\right\|_F=\|\boldsymbol{A}\|_F \text { and }\left\|\boldsymbol{A}^*\right\|_2=\|\boldsymbol{A}\|_2 \text {. }
```
"""

# ╔═╡ 4cea17b5-b132-4de2-9294-db3ebeccbbd3
cm"""
$(define("Absolute Norm"))
A vector norm on ``\mathbb{C}^n`` is an absolute norm if ``\|\boldsymbol{x}\|=\||\boldsymbol{x}|\|`` for all ``\boldsymbol{x} \in \mathbb{C}^n``. Here ``|\boldsymbol{x}|:=\left[\left|x_1\right|, \ldots,\left|x_n\right|\right]^T``, the absolute values of the components of ``\boldsymbol{x}``. 
$(ebl())

$(bbl("Remarks",""))
- The vector ``p`` norms are absolute norms. 
- A vector norm on ``\mathbb{C}^n`` is an absolute norm if and only if it is a monotone norm, i.e.,
```math
\left|x_i\right| \leq\left|y_i\right|, i=1, \ldots, n \Longrightarrow\|\boldsymbol{x}\| \leq\|\boldsymbol{y}\| \text {, for all } \boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^n \text {. }
```
"""

# ╔═╡ 34328fe5-caf3-441c-8957-0190ecc62170
cm"""
$(bth("8.7 (Perturbation in the Right-Hand Side)"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is nonsingular, ``\boldsymbol{b}, \boldsymbol{e} \in \mathbb{C}^n, \boldsymbol{b} \neq \mathbf{0}`` and ``\boldsymbol{A x}=\boldsymbol{b}, \boldsymbol{A} \boldsymbol{y}=\boldsymbol{b}+\boldsymbol{e}``. Then
```math
\frac{1}{K(\boldsymbol{A})} \frac{\|\boldsymbol{e}\|}{\|\boldsymbol{b}\|} \leq \frac{\|\boldsymbol{y}-\boldsymbol{x}\|}{\|\boldsymbol{x}\|} \leq K(\boldsymbol{A}) \frac{\|\boldsymbol{e}\|}{\|\boldsymbol{b}\|},
```
where ``K(\boldsymbol{A})=\|\boldsymbol{A}\|\left\|\boldsymbol{A}^{-1}\right\|`` is the __condition number of ``\boldsymbol{A}``__.
"""

# ╔═╡ 17f6efeb-ccbb-4aa0-bf02-056195f91604
cm"""
$(bth("8.8 (Perturbation and Residual)"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}, \boldsymbol{b} \in \mathbb{C}^n``, ``\boldsymbol{A}`` is nonsingular and ``\boldsymbol{b} \neq \mathbf{0}``. Let ``\boldsymbol{r}(\boldsymbol{y})=\boldsymbol{A} \boldsymbol{y}-\boldsymbol{b}`` for ``\boldsymbol{y} \in \mathbb{C}^n``. If ``\boldsymbol{A} \boldsymbol{x}=\boldsymbol{b}`` then
```math
\frac{1}{K(\boldsymbol{A})} \frac{\|\boldsymbol{r}(\boldsymbol{y})\|}{\|\boldsymbol{b}\|} \leq \frac{\|\boldsymbol{y}-\boldsymbol{x}\|}{\|\boldsymbol{x}\|} \leq K(\boldsymbol{A}) \frac{\|\boldsymbol{r}(\boldsymbol{y})\|}{\|\boldsymbol{b}\|} .
```
"""

# ╔═╡ 9ee3d314-acee-4646-8139-6fd1706e7292
cm"""
$(bth("8.9 (Perturbation in Matrix)"))
Suppose ``\boldsymbol{A}, \boldsymbol{E} \in \mathbb{C}^{n \times n}, \boldsymbol{b} \in \mathbb{C}^n`` with ``\boldsymbol{A}`` nonsingular and ``\boldsymbol{b} \neq \mathbf{0}``. If ``r:=\left\|\boldsymbol{A}^{-1} \boldsymbol{E}\right\|<1`` then ``\boldsymbol{A}+\boldsymbol{E}`` is nonsingular. If ``\boldsymbol{A} \boldsymbol{x}=\boldsymbol{b}`` and ``(\boldsymbol{A}+\boldsymbol{E}) \boldsymbol{y}=\boldsymbol{b}`` then
```math
\begin{aligned}
& \frac{\|\boldsymbol{y}-\boldsymbol{x}\|}{\|\boldsymbol{y}\|} \leq r \leq K(\boldsymbol{A}) \frac{\|\boldsymbol{E}\|}{\|\boldsymbol{A}\|} \\
& \frac{\|\boldsymbol{y}-\boldsymbol{x}\|}{\|\boldsymbol{x}\|} \leq \frac{r}{1-r} \leq \frac{K(\boldsymbol{A})}{1-r} \frac{\|\boldsymbol{E}\|}{\|\boldsymbol{A}\|}
\end{aligned}
```
"""

# ╔═╡ de24c4e8-65c8-4952-a87f-5196ef74891a
cm"""
$(bth("8.10 (Spectral Condition Number)"))
Suppose ``\boldsymbol{A} \in \mathbb{C}^{n \times n}`` is nonsingular with singular values ``\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n>0`` and eigenvalues ``\left|\lambda_1\right| \geq\left|\lambda_2\right| \geq \cdots \geq`` ``\left|\lambda_n\right|>0``. Then
```math
K_2(\boldsymbol{A})= \begin{cases}\lambda_1 / \lambda_n, & \text { if } \boldsymbol{A} \text { is positive definite }, \\ \left|\lambda_1\right| /\left|\lambda_n\right|, & \text { if } \boldsymbol{A} \text { is normal }, \\ \sigma_1 / \sigma_n, & \text { in general. }\end{cases}
```
"""

# ╔═╡ db507c85-6847-4e29-aaf5-7314b22abe2d
cm"""
$(define("9.1 (Least Squares Problem (LSQ))"))
Suppose ``m, n \in \mathbb{N}, \boldsymbol{A} \in \mathbb{C}^{m \times n}`` and ``\boldsymbol{b} \in \mathbb{C}^m``. To find ``\boldsymbol{x} \in \mathbb{C}^n`` that minimizes ``E: \mathbb{C}^n \rightarrow \mathbb{R}`` given by
```math
E(\boldsymbol{x}):=\|\boldsymbol{A} \boldsymbol{x}-\boldsymbol{b}\|_2^2,
```
is called the __least squares problem__. A minimizer ``\boldsymbol{x}`` is called a __least squares solution__.
"""

# ╔═╡ 083fd2fa-582a-414f-87b3-012a3509a049
cm"""
We show the following results, valid for any ``m, n \in \mathbb{N}, \boldsymbol{A} \in \mathbb{C}^{m \times n}`` and ``\boldsymbol{b} \in \mathbb{C}^n``.

$(bth("9.1 (Existence)")) The least squares problem always has a solution.
$(ebl())

$(bth("9.2 (Uniqueness)")) The solution of the least squares problem is unique if and only if ``\boldsymbol{A}`` has linearly independent columns.
$(ebl())

$(bth("9.3 (Characterization)"))
``\boldsymbol{x} \in \mathbb{C}^n`` is a solution of the least squares problem if and only if ``\boldsymbol{A}^* \boldsymbol{A x}=\boldsymbol{A}^* \boldsymbol{b}``.
$(ebl())
"""

# ╔═╡ 71ea919d-fcf1-4282-8d71-2888406d31d0
cm"""
$(bbl("Example 9.2", "(Linear Regression)"))
We want to fit a straight line ``p(t)=x_1+x_2 t`` to ``m \geq 2`` given data ``\left(t_k, y_k\right) \in \mathbb{R}^2, k=1, \ldots, m``. This is part of the linear regression process in statistics. We obtain the linear system
```math
\boldsymbol{A} \boldsymbol{x}=\left[\begin{array}{c}
p\left(t_1\right) \\
\vdots \\
p\left(t_m\right)
\end{array}\right]=\left[\begin{array}{ll}
1 & t_1 \\
\vdots \\
1 & t_m
\end{array}\right]\left[\begin{array}{l}
x_1 \\
x_2
\end{array}\right]=\left[\begin{array}{c}
y_1 \\
\vdots \\
y_m
\end{array}\right]=\boldsymbol{b}
```
"""

# ╔═╡ bdce2e27-5fe4-47cc-b27e-7cfc2e320a45
cm"""
$(ex())
Give the data,

| ``t`` | ``1.0`` | ``2.0`` | ``3.0`` | ``4.0`` |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| y | 3.1 |  1.8 | 1.0 |  0.1 |

Find the __least squares line__ (Linear Regression)
"""

# ╔═╡ 4e890dc7-db6a-4213-ac44-a2bb80e542c6
cm"""
$(bth("(Spectral Condition Number)"))
of ``\boldsymbol{A}^* \boldsymbol{A}`` 

Suppose ``1 \leq n \leq m`` and that ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` has linearly independent columns. Then
```math
K_2\left(\boldsymbol{A}^* \boldsymbol{A}\right):=\left\|\boldsymbol{A}^* \boldsymbol{A}\right\|_2\left\|\left(\boldsymbol{A}^* \boldsymbol{A}\right)^{-1}\right\|_2=\frac{\lambda_1}{\lambda_n}=\frac{\sigma_1^2}{\sigma_n^2}=K_2(\boldsymbol{A})^2,
```
where ``\lambda_1 \geq \cdots \geq \lambda_n>0`` are the eigenvalues of ``\boldsymbol{A}^* \boldsymbol{A}``, and ``\sigma_1 \geq \cdots \geq \sigma_n>0`` are the singular values of ``\boldsymbol{A}``.

Proof Since ``\boldsymbol{A}^* \boldsymbol{A}`` is Hermitian it follows from Theorem 8.10 that ``K_2(\boldsymbol{A})=\frac{\sigma_1}{\sigma_n}`` and ``K_2\left(\boldsymbol{A}^* \boldsymbol{A}\right)=\frac{\lambda_1}{\lambda_n}``. But ``\lambda_i=\sigma_i^2`` by Theorem 7.2 and the proof is complete.
"""

# ╔═╡ 7ad74395-20ad-4d02-92a9-788790382fa4
cm"""
$(ex())
Solve ``Ax=b`` where
```math
\boldsymbol{A}=\left[\begin{array}{rrr}1 & 3 & 1 \\ 1 & 3 & 7 \\ 1 & -1 & -4 \\ 1 & -1 & 2\end{array}\right] \quad \text{and}\quad \boldsymbol{b}=\left[\begin{array}{l}1 \\ 1 \\ 1 \\ 1\end{array}\right]
```

"""

# ╔═╡ 1b72137a-1d39-49bc-9d2f-8d60262dccb1
cm"""
$(bth("(The Generalized Inverse)"))
For any ``m, n \in \mathbb{N}`` and any ``\boldsymbol{A} \in \mathbb{C}^{m \times n}`` there is a unique matrix ``\boldsymbol{A}^{\dagger} \in \mathbb{C}^{n \times m}`` such that
```math
\boldsymbol{A} \boldsymbol{A}^{\dagger} \boldsymbol{A}=\boldsymbol{A}, \boldsymbol{A}^{\dagger} \boldsymbol{A} \boldsymbol{A}^{\dagger}=\boldsymbol{A}^{\dagger},\left(\boldsymbol{A}^{\dagger} \boldsymbol{A}\right)^*=\boldsymbol{A}^{\dagger} \boldsymbol{A},\left(\boldsymbol{A} \boldsymbol{A}^{\dagger}\right)^*=\boldsymbol{A} \boldsymbol{A}^{\dagger}\tag{👽}
```

If ``\boldsymbol{U}_1 \boldsymbol{\Sigma}_1 \boldsymbol{V}_1^*`` is a singular value factorization of ``\boldsymbol{A}`` then
```math
\boldsymbol{A}^{\dagger}=\boldsymbol{V}_1 \boldsymbol{\Sigma}_1^{-1} \boldsymbol{U}_1^*
```
"""

# ╔═╡ 2a6ad2c4-16f2-4a3c-9bf8-b198fe560626
cm"""
$(bth("(Orthogonal Projections)"))
Given ``m, n \in \mathbb{N}, \boldsymbol{A} \in \mathbb{C}^{m \times n}`` of rank ``r``, and let ``\mathcal{S}`` be one of the subspaces ``\mathcal{R}(\boldsymbol{A}), \mathcal{N}\left(\boldsymbol{A}^*\right)``. The orthogonal projection of ``\boldsymbol{v} \in \mathbb{C}^m`` into ``\mathcal{S}`` can be written as a matrix ``\boldsymbol{P}_{\mathcal{S}}`` times the vector ``\boldsymbol{v}`` in the form ``\boldsymbol{P}_{\mathcal{S}} \boldsymbol{v}``, where
```math
\begin{gathered}
\boldsymbol{P}_{\mathcal{R}(\boldsymbol{A})}=\boldsymbol{A} \boldsymbol{A}^{\dagger}=\boldsymbol{U}_1 \boldsymbol{U}_1^*=\sum_{j=1}^r \boldsymbol{u}_j \boldsymbol{u}_j^* \in \mathbb{C}^{m \times m} \\
\boldsymbol{P}_{\mathcal{N}\left(\boldsymbol{A}^*\right)}=\boldsymbol{I}-\boldsymbol{A} \boldsymbol{A}^{\dagger}=\boldsymbol{U}_2 \boldsymbol{U}_2^*=\sum_{j=r+1}^m \boldsymbol{u}_j \boldsymbol{u}_j^* \in \mathbb{C}^{m \times m} .
\end{gathered}
```
where ``\boldsymbol{A}^{\dagger}`` is the generalized inverse of ``\boldsymbol{A}``, and ``\boldsymbol{A}=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^* \in \mathbb{C}^{m \times n}`` is a singular value decomposition of ``\boldsymbol{A}(c f``. (9.7)).
"""

# ╔═╡ 16ff57d8-3c79-4442-80f0-f8a43a949649
cm"""
$(bbl("Collary","(LSQ Characterization Using Generalized Inverse)")) ``\boldsymbol{x} \in \mathbb{C}^n`` solves the least squares problem ``\min _x\|\boldsymbol{A x}-b\|_2^2`` if and only if ``\boldsymbol{x}=\boldsymbol{A}^{\dagger} \boldsymbol{b}+\boldsymbol{z}``, where ``\boldsymbol{A}^{\dagger}`` is the generalized inverse of ``\boldsymbol{A}`` and ``\boldsymbol{z} \in \mathcal{N}(\boldsymbol{A})``.
"""

# ╔═╡ 04a47b4d-9a14-4225-a907-a63ac6ee70d6
cm"""
$(bth("(Minimal Norm Solution)"))
The least squares solution with minimal Euclidian norm is ``\boldsymbol{x}=\boldsymbol{A}^{\dagger} \boldsymbol{b}`` corresponding to ``\boldsymbol{z}=\mathbf{0}``.
"""

# ╔═╡ b382223f-9c78-4671-97ec-d898654dfd0a
cm"""
$(bth("(Splitting Matrices for J, GS and SOR)")) The splitting matrices ``\boldsymbol{M}_J, \boldsymbol{M}_1`` and ``\boldsymbol{M}_\omega`` for the J, GS and SOR methods are given by
```math
\boldsymbol{M}_J=\boldsymbol{D}, \quad \boldsymbol{M}_1=\boldsymbol{D}-\boldsymbol{A}_L, \quad \boldsymbol{M}_\omega=\omega^{-1} \boldsymbol{D}-\boldsymbol{A}_L
```
"""

# ╔═╡ c03160d4-b063-4899-8df6-6ebd66cc4617
cm"""
$(ex())
Solve 
```math
\begin{array}{llllllllll} 
& 10 x_1&-&x_2&+&2 x_3&& & =&6 \\ 
-&x_1&+&11 x_2&-&x_3&+&3 x_4 & =&25 \\ 
&2 x_1&-&x_2&+&10 x_3&-&x_4 & =&-11 \\ 
&&&3 x_2&-&x_3&+&8 x_4 & =&15
\end{array}
```
"""

# ╔═╡ 658ced39-70ff-4fb1-8fdd-b314aca309ce
cm"""
In the following we assume that ``\boldsymbol{G} \in \mathbb{C}^{n \times n}, \boldsymbol{c} \in \mathbb{C}^n`` and ``\boldsymbol{I}-\boldsymbol{G}`` is nonsingular. We let ``\boldsymbol{x} \in \mathbb{C}^n`` be the unique fixed point satisfying ``\boldsymbol{x}=\boldsymbol{G} \boldsymbol{x}+\boldsymbol{c}``.

$(define("(Convergence)"))
We say that the iterative method ``\boldsymbol{x}_{k+1}:=\boldsymbol{G} \boldsymbol{x}_k+\boldsymbol{c}`` converges if the sequence ``\left\{\boldsymbol{x}_k\right\}`` converges for any starting vector ``\boldsymbol{x}_0``.
$(ebl())

$(bth("(Convergence of an Iterative Method)"))
The iterative method ``\boldsymbol{x}_{k+1}:=\boldsymbol{G} \boldsymbol{x}_k+\boldsymbol{c}`` converges if and only if ``\lim _{k \rightarrow \infty} \boldsymbol{G}^k=\mathbf{0}``.
$(ebl())

$(bth("(Sufficient Condition for Convergence)"))
If ``\|\boldsymbol{G}\|<1`` then the iteration ``\boldsymbol{x}_{k+1}=\boldsymbol{G} \boldsymbol{x}_k+\boldsymbol{c}`` converges.
$(ebl())
"""

# ╔═╡ 8c4bf470-1408-4057-8d7e-216b1b38c1fa
cm"""
$(bth("(When Does an Iterative Method Converge?)")) Suppose ``\boldsymbol{G} \in \mathbb{C}^{n \times n}`` with ``\boldsymbol{I}-\boldsymbol{G}`` nonsingular and let ``\boldsymbol{c} \in \mathbb{C}^n``. The iteration ``\boldsymbol{x}_{k+1}=\boldsymbol{G} \boldsymbol{x}_k+\boldsymbol{c}`` converges if and only if ``\rho(\boldsymbol{G})<1`` where
```math
\rho(\boldsymbol{G}):=\max _{\lambda \in \sigma(\boldsymbol{G})}|\lambda|.
```
and is called __spectral radius__ of ``G``.
"""

# ╔═╡ 6c86adaf-2850-46f1-b421-e4ebc97e6631
cm"""
$(bbl("Lemma","(Be Careful When Stopping)")) If ``\boldsymbol{x}_k=\boldsymbol{G} \boldsymbol{x}_{k-1}+\boldsymbol{c}, \boldsymbol{x}=\boldsymbol{G} \boldsymbol{x}+\boldsymbol{c}`` and ``\|\boldsymbol{G}\|<1`` then
```math
\left\|\boldsymbol{x}_k-\boldsymbol{x}_{k-1}\right\| \geq \frac{1-\|\boldsymbol{G}\|}{\|\boldsymbol{G}\|}\left\|\boldsymbol{x}_k-\boldsymbol{x}\right\|, \quad k \geq 1
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QRCoders = "f42e9828-16f3-11ed-2883-9126170b272d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Colors = "~0.12.11"
Combinatorics = "~1.0.2"
CommonMark = "~0.8.12"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Latexify = "~0.16.5"
PlotThemes = "~3.2.0"
Plots = "~1.40.5"
PlutoExtras = "~0.7.12"
PlutoUI = "~0.7.59"
PrettyTables = "~2.3.2"
QRCoders = "~1.4.5"
Roots = "~2.1.6"
SymbolicUtils = "~2.1.2"
Symbolics = "~5.35.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "057e33cad75fbc8a952cd7fbb187792e23e4affc"

[[deps.ADTypes]]
git-tree-sha1 = "eea5d80188827b35333801ef97a40c2ed653b081"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.9.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "PrecompileTools"]
git-tree-sha1 = "3faae67b8899797592335832fccf4b3c80bb04fa"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.15"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "d7477ecdafb813ddee2ae727afa94e9dcb5f3fb0"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.112"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "490392af2c7d63183bfa2c8aaa6ab981c5ba7561"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.14"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "30a1848c4f4fc35d1d4bbbd125650f6a11b5bc6c"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.7"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.Extents]]
git-tree-sha1 = "81023caa0021a41712685887db1fc03db26f41f5"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.4"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "62ca0547a14c57e98154423419d8a342dca75ca9"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.4"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "59107c179a586f0fe667024c5eb7033e81333271"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.2"

[[deps.GeoInterface]]
deps = ["Extents", "GeoFormatTypes"]
git-tree-sha1 = "2f6fce56cdb8373637a6614e14a5768a88450de2"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.7"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d65554bad8b16d9562050c67e7223abf91eaba2f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.13+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "e459fda6b68ea8684b3fcd513d2fd1e5130c4402"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.16.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "fa7fd067dca76cadd880f1ca937b4f387975a9f5"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.16.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.MarchingCubes]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "301345b808264ae42e60d10a519e55c5d992969b"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3eba928678787843e504c153a9b8e80d7d73ab17"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.5.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "f4cb457ffac5f5cf695699f82c537073958a6a6c"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.2+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "DocStringExtensions", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoUI", "REPL", "Random"]
git-tree-sha1 = "681f89bdd5c1da76b31a524af798efb5eb332ee9"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.13"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QRCoders]]
deps = ["FileIO", "ImageCore", "ImageIO", "ImageMagick", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "b3e5fcc7a7ade2d43f0ffd178c299b7a264c268a"
uuid = "f42e9828-16f3-11ed-2883-9126170b272d"
version = "1.4.5"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "b034171b93aebc81b3e1890a036d13a9c4a9e3e0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "ChainRulesCore", "CommonSolve", "Printf"]
git-tree-sha1 = "48a7925c1d971b03bb81183b99d82c1dc7a3562f"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.1.8"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "98ca7c29edd6fc79cd74c61accb7010a4e7aee33"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.6.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "50ed64cd5ad79b0bef71fdb6a11d10c3448bfef0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.56.1"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "e39c5f217f9aca640c8e27ab21acf557a3967db5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.10"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"
weakdeps = ["Adapt", "GPUArraysCore", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "0225f7c62f5f78db35aae6abb2e5cabe38ce578f"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.31"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fb099adbd7504f1e68b4512828e9d94197a8b889"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "a7ecd8fb3ee77a47ac033395069117763469c49e"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "2.1.3"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "ForwardDiff", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "653f8337e6605b4a84e5da15add78f1de8f4bae5"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.35.1"

    [deps.Symbolics.extensions]
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxCoreExt = "LuxCore"
    SymbolicsPreallocationToolsExt = "PreallocationTools"
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LuxCore = "bb33d45b-7691-41d6-9220-0943567d0623"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "6f0cee95e74d1f6891ba6b35b8b219fd3d11b567"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.4.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "38f139cc4abf345dd4f22286ec000728d5e8e097"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3a6f063d690135f5c1ba351412c82bae4d1402bf"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.25"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorTypes", "Contour", "Crayons", "Dates", "FileIO", "FreeTypeAbstraction", "LazyModules", "LinearAlgebra", "MarchingCubes", "NaNMath", "Printf", "SparseArrays", "StaticArrays", "StatsBase", "Unitful"]
git-tree-sha1 = "ae67ab0505b9453655f7d5ea65183a1cd1b3cfa0"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "2.12.4"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "7dfa0fd9c783d3d0cc43ea1af53d69ba45c447df"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─06790633-b8aa-45ba-ac67-c416c88b166a
# ╟─c6241f8a-60ff-44c0-b363-f2d91b2c5eb0
# ╟─01183bde-075c-4848-995a-b23ffeeb97c8
# ╟─a8e556ea-48ca-4a43-914a-7171a6491186
# ╟─73b56b54-22a2-4fa0-8eed-ab8a23cebc74
# ╟─06de9303-589d-40cf-ac95-4df2020af3a6
# ╠═c59dffb4-c419-461b-8096-e27171be0a87
# ╟─a8948e17-2846-431f-9765-6359eaeb20a9
# ╟─db84f278-2b61-40fd-b0a8-bb132cff5f18
# ╟─df6c46d3-fa32-4709-a392-463167b33c46
# ╠═b976040e-4974-4bfc-a6a9-e3926a1f2eef
# ╠═5356f40b-2cc7-4490-9752-115cae126839
# ╠═e3e2f379-f7bb-4ccf-85e0-378ff479f3de
# ╠═6ccfccef-09cc-42d1-a5c8-3899801b4438
# ╟─6dec69df-0c80-4790-b5c5-12fbe2dc41b8
# ╠═17263bb2-6964-48a6-b31e-b20f168f35a2
# ╠═30bc89ef-58af-4bed-b172-aa6b4e2c8491
# ╟─7e88c07e-394e-46eb-b1e5-89fb65287e36
# ╟─21a9d1b8-fdb4-44e6-99da-f97ac172a9a3
# ╟─5742ff61-39a1-41e6-9136-1f148516d5e0
# ╟─c7f1e8cd-da40-469f-b8cf-42869d5b46ed
# ╟─7f18dc1e-e040-4f0d-b3bf-a5477489a1ab
# ╟─ffcc0c43-bbb5-433c-a98f-b8421a846485
# ╟─8aa563a5-8264-4bf0-89c9-c1daa74bd4d6
# ╠═dbf3b4de-06ce-46ff-9178-cc28961ab3e5
# ╠═5ab2ce16-f52d-4107-9035-4dc47df19fcd
# ╟─0fb73bc3-e70d-4de2-a8fa-49938c6f3a60
# ╠═4dbfa5bd-dfcc-4195-8337-02f8bed8748a
# ╠═49ee7daf-8b09-4303-8aa6-05ae57977688
# ╟─21e63aa1-c0de-4980-b973-9afd6b1dde87
# ╟─65875093-4ff0-4684-a355-678d727f36f7
# ╟─6d79d33f-8caa-4fa9-b1e3-baf6bb25d519
# ╟─374034ce-15d5-403d-9d2a-aa760cc0a269
# ╟─e06e478f-b2a9-4b33-8a36-3bc87e53f62a
# ╟─8d801c77-08e3-45c0-b27c-0d36f0bfb29e
# ╟─0beeed3b-905e-4ac1-b19c-a06cec3f160f
# ╟─41e2ac3f-5488-4506-9278-b7857d8ad4e7
# ╟─d5e9beb5-56a9-4a99-903d-719dc18cf710
# ╟─7513851f-48f8-4e2f-8c97-baee10e9e4c7
# ╟─9a95f3c9-0e47-422a-9b0e-35fc8c528308
# ╟─918ca6c0-ac2d-4072-b7d9-b73de95248fd
# ╠═a5240ee4-1f4d-4229-95f4-fe115645c92c
# ╠═05dfaf81-8b67-43de-b249-f1e23c3664f2
# ╟─77b036e6-7ddb-4a8e-8da5-3e0d1c9305ca
# ╟─3dcc9b0b-2daa-4695-9e77-311672ac511b
# ╟─1a8d7932-4363-4749-ab31-a8ecf5288169
# ╟─4486fd72-e3f4-4e78-bc57-db8474dca974
# ╟─50ac8333-8f68-4f85-be14-dc99ba7fc398
# ╟─3f43de8a-b6b4-4d53-aad4-bd3356f2d2a8
# ╟─519206e2-863d-4099-9f20-11cdad64a5e0
# ╟─29f4fc5f-9f93-4004-9089-061fc14517f9
# ╟─1740ab31-1c55-4bfc-a34e-dc10852de7ff
# ╟─45ef1c9d-1052-4656-b4e5-1704c39013ee
# ╟─cda6aea6-61ca-483f-98f9-2274ef50e049
# ╟─0ec8c4e5-ba42-48c9-95ef-47457451e372
# ╟─adc9275a-493c-40c8-9b9a-05d38e43036e
# ╟─58451143-9990-4426-9d42-50a5d3d9849b
# ╟─5fef6c3a-9809-4906-b5b8-4d00c72439a8
# ╟─da8fd297-56f4-4175-8694-093de356dff9
# ╟─71a38e11-f634-4e57-a277-e9eecd324381
# ╠═8daf8b4c-39ed-47d2-932d-28b4c31ad6f7
# ╟─b5c0a76a-dc08-46f9-9b2d-f1b528f1aa87
# ╠═420504bd-7c7b-4a09-a0c1-203b36ed340e
# ╠═e9ab9e3f-47ae-48e7-ae3d-4afe1618919a
# ╟─e141f527-3783-4293-a890-10c05aa166dd
# ╟─6d5bf021-a49e-4726-a446-900422cc6703
# ╠═6ed68d7b-7c50-4150-86f0-905fe6fd6a25
# ╟─a4226ef9-8d8f-4508-9e88-6dc09f6410be
# ╟─dff54fcc-2793-41d9-826a-8321e536b7d7
# ╟─953f4449-9dff-4275-b32b-075d727d2d2e
# ╟─2e68ed26-ae75-42f7-bc5d-fb23348e6711
# ╟─6897cf15-2370-4114-a574-ea137096a175
# ╟─d312c0d7-f80b-4da7-9ae3-9cdb2f783635
# ╠═9521e54d-5a82-44f8-a4b6-94b229a58ec5
# ╟─a0560ea2-1ab5-425a-a418-92e196b51f92
# ╟─b9a41148-048e-463c-b0c1-85de79ec6ee1
# ╟─aa5e1849-e9e4-4b2f-90b7-9716ad441a85
# ╟─82c571b9-23a6-48b6-b58c-a42269bfa429
# ╟─8816f81b-3f55-49e0-ae40-ff6672c951a9
# ╟─2d53670e-772e-41b3-9076-3193b40a5c09
# ╟─12e05a32-7a96-420c-a90b-1b76748cb09f
# ╟─d068efed-78d5-4913-a96c-3ec28e5988f7
# ╟─b96f87dd-14b9-49ce-970a-b5093fc13e94
# ╟─e0504a1a-d598-4ca7-9023-77d5d57814e0
# ╠═a78b0614-3a52-4fb8-acc7-147707e1d456
# ╟─fd98de8b-bb26-4c7c-8b65-d16f3fe57a09
# ╟─7f2f4c83-63c9-401f-b6c5-25722674a68c
# ╟─c4ec32de-d6e1-4ef4-b870-cd5c5e7f77fa
# ╟─7090fc3e-97e5-4819-af88-1bb2dd51d7ec
# ╟─802688ba-a8f5-419b-a135-b6525679fb42
# ╠═d5dd27d9-16dc-44b0-9128-63e16cc9db10
# ╟─93c8a52c-e3d9-42a5-b9be-11e6a090fa55
# ╟─c9344f33-0372-48de-ae18-819e05d5852f
# ╟─81cb8f6c-d82b-41eb-bd06-f4a601954785
# ╟─91739d93-dfda-4599-9c4b-68c6ee7f8df1
# ╟─e440727e-955e-4db6-a9d5-c1125917f56e
# ╟─fb35d4ab-fc81-4618-be09-4a2911ff4566
# ╟─2ad11900-f756-4600-9565-55108bd0e296
# ╟─3b933d21-8ea9-4790-8588-4662a84481ba
# ╟─43ffc137-855c-4e16-8ebd-ce2f73ad7161
# ╟─e90b24e7-1c5a-4190-9af6-34fe939ad47c
# ╠═aee4e1e3-7ce8-48a9-b258-b0ab46d52925
# ╟─32228ef3-45bb-473c-9a04-21070b475b19
# ╟─380db7d1-9bb6-44e3-9063-6d98bceaf15f
# ╟─f7f904ea-2c8b-473f-a37f-53e6ed627a77
# ╟─5225f24c-d76d-4011-8fd3-8c9e16e886c6
# ╟─faf07e34-58fd-4be3-a81d-9fddd668e582
# ╟─ea4df632-2f6f-4cde-acf7-348b8ebbfb63
# ╟─3a9f7288-bcac-451e-b361-86102a018c3f
# ╟─dc85a1f4-bff2-466b-98a8-44328d6cf7c1
# ╟─e9ea8099-acfd-4df4-a19a-f510888e1f98
# ╟─0275d178-b918-445a-be9e-eac9c7aaa471
# ╟─f68a295a-8fe6-47d1-9ca6-9b6199163ed7
# ╟─6d988a2c-f6ae-418d-9cd4-11207ec7a2ce
# ╠═02c7ec2d-e57f-49d5-a8c3-a6f0637e026f
# ╟─27e009e2-2e77-4e58-a2d8-50a9efb99e32
# ╟─9579da39-7241-4266-918a-0c46b656cbf3
# ╟─6e8ecdb1-7d32-429f-9d23-8b79db0cc75d
# ╟─41526c5c-0be0-4bef-af4c-557c28b918cd
# ╟─694c5383-28cf-4e69-90bc-366a40c1e0a9
# ╟─a3a98e03-e9ea-424e-972d-f8367cd52642
# ╟─35e8b2a4-06d7-4fa4-9b02-63cdffaf1961
# ╟─c0f9964e-b8f8-4eca-8126-64a139e20afc
# ╟─b66c037a-d1cc-42d2-a45c-9f0190b8f28d
# ╟─ff626183-36f2-455a-b141-1b9725219f41
# ╟─b012a15b-be75-4cb4-96ff-d679ed4376d1
# ╟─e42b6939-e52d-4b9b-adba-02f4549fd08b
# ╟─e53a67cc-c38f-4006-8d20-7bbc17e26f66
# ╠═8509b4b9-b4c8-47e5-84ae-c21b47eafa16
# ╟─57a9ebc0-a35f-4900-ba3f-9a7df19ecfbf
# ╟─6de81e3d-53b3-458c-a9cf-edcef35e0db3
# ╟─2bebb2e6-55b6-4a33-9409-2b901b6dfc5e
# ╟─5e8ea837-cd68-4dd1-9d4b-925223fc63f9
# ╟─a4a31455-1959-4af0-9b56-ec9dec4c94d5
# ╟─489d96e6-14b8-4096-829f-71059ad6d25c
# ╟─614c7edf-1133-4e93-a5d0-77747b840ca7
# ╟─17312234-ed45-4e26-868e-3f25c14f73bd
# ╟─310f6fa8-f185-4901-8368-b2b268e40bca
# ╟─ef18b323-485f-48f4-95d0-9063bb6ef2e1
# ╠═3d678fbd-65f7-4516-a684-6724c66970de
# ╟─a916bede-5304-4120-a929-979d8fbff63c
# ╠═b5ba03bb-8c78-4a14-91d8-f4b110db8885
# ╟─5ae18425-3dde-487d-80f7-127d59a18bbb
# ╟─e28d00ed-710d-4c08-af2e-0f9562d64be2
# ╟─4eade3b9-589f-492a-9965-03eb74afd493
# ╟─7627739e-9f51-4197-8ece-ad3e17b0f906
# ╟─a6148859-78c2-4c40-aa0e-0373870a74b8
# ╟─c5077118-527a-448b-bfbc-ad5f00082b7b
# ╟─9d950ca7-0162-41a6-b59f-ba2d64d48f4e
# ╠═7e6ce1a2-8de1-4170-9585-62b25311e646
# ╠═d94ca4da-f011-4c8c-b92e-985d29c4f3e5
# ╟─0865697a-e057-4be3-9aa9-a1bf8831829d
# ╠═7c3661da-0ee8-4877-8d50-07c6704dbbb0
# ╟─8acef7cd-a100-46be-b3ce-6b00ab56f479
# ╟─a71dc775-4826-494d-8b83-62274561e6be
# ╟─864a50be-eabc-44be-891e-64d7caa2d8ef
# ╟─7fd8fae9-c083-46c7-a0b9-1c8317641669
# ╟─b35a0a84-880b-4da3-a01d-f7d7fd1bb2a4
# ╟─2d60f5a5-139f-4451-af2e-33ee04355709
# ╠═dca4ccc0-bb15-41b3-9fb5-a3d74bf552c9
# ╠═f79e868a-5ed7-4439-ac3c-9229f64360d2
# ╟─04b6d34a-7a30-4337-a862-a84f17680ac6
# ╠═3ba9e235-869c-4011-afeb-c09385260140
# ╠═3e8a3007-5161-40f5-bf7d-47d8aeecb1aa
# ╟─7ea5927a-81bb-404a-aa62-67bce4e00313
# ╟─211f7f09-31f4-44d5-8e8f-6bb5a2920ae7
# ╟─188ddb61-6cba-4485-83a0-ff37870cebed
# ╟─40820595-a553-4e90-9b78-d6d4b3c473ae
# ╠═ca414318-5a80-423e-b4a3-561e282e4111
# ╟─37b9d55e-bc0d-438e-bb99-690b056bd2df
# ╟─74d2a314-2eae-491e-a4ad-dc1184195c00
# ╟─be674a64-cb53-47d0-bdfb-33d72266335d
# ╠═2bb49eb8-1c93-44f3-a05b-745594830356
# ╟─cc72a9e5-1f40-496e-8f87-79cf6734da89
# ╟─2e9130ea-3634-45a3-9658-d7b160ed491d
# ╠═fbf9565f-f1d6-4a50-b81d-4863eff22d3b
# ╠═c6576e1e-f86f-46a4-bd2f-7b82550c231d
# ╟─22cac4d7-7c2a-43e9-8d94-ed6df2508fed
# ╟─f0a58c35-47dc-4666-befc-08c502e6e229
# ╟─5a573a95-6092-4def-b900-1b80c3aca31f
# ╟─3f925be3-06da-4b91-b9c8-1749ba55b3d0
# ╟─37344680-8fd2-4440-9f03-b7b325c8965c
# ╟─ea67c00c-4aed-46ff-a11f-5455553901c0
# ╟─e51e2e47-e389-4afa-b1e0-1ed5bf3768f6
# ╟─bfaf9649-e08a-4adf-9794-539277902565
# ╟─e7a24ffb-1a85-433c-a9c2-4868f9ba4205
# ╟─458a0966-e841-48dc-8e04-6f617fb70b6b
# ╟─2aa4b07a-f768-4986-814d-eb5c6a279adc
# ╟─bf91a78d-429c-46c1-8a09-8f7551c32d6a
# ╟─823551c5-01bc-4217-82d9-a87bbb7f03fd
# ╟─4f2e0531-99fa-4c48-867b-710c80b42be5
# ╟─9fdb84a6-2576-4b58-bf2e-6daa75434f1f
# ╟─e4782556-d145-4b88-ba80-359da77b2362
# ╟─d93f5977-87d5-42cf-b383-8c861a69ccf5
# ╠═1660650f-917c-4354-acb7-c287f53548b8
# ╠═8221aa32-82ae-4ec6-b40b-0ec7c54f92de
# ╟─3ea9b542-6182-48d3-9fc6-29b7d1c28920
# ╟─2c711e97-7f1f-4c7c-b7c0-3cfb56e912fd
# ╟─a657c07f-ab42-4051-8e03-fedb8af62f6b
# ╟─b5c411da-2629-4942-944b-5f3f4245b74e
# ╟─43dc89bb-3c51-4aea-98ee-4df58ccd2997
# ╟─db3892ca-816a-4982-adb0-ef09109a99a0
# ╠═df706f98-7cf7-4b84-9868-0e67e9e4c34a
# ╟─ce0f671e-395d-4cd4-a95b-5ec67e7c1fd0
# ╠═6e123b71-f3bd-40b5-ad65-82515b0a44df
# ╟─63c34957-6379-4723-b537-f9753f6fe814
# ╟─4830028b-2e80-48b1-8dc0-79597852b828
# ╟─6be7cbc4-023b-4d2e-bda5-36a996fad82c
# ╟─2f9a1b80-5247-4af9-bf89-15a50eb7e5c6
# ╟─70686047-be6a-4d77-8cc3-dcb38b7bc378
# ╟─afa0532c-de60-44d4-a166-befcdc1959e4
# ╟─b23fc699-d1e8-4c83-8e61-06dc4ab8ea6a
# ╟─a2dc9322-049b-4437-a29e-1cd23f816059
# ╟─4bf58173-2848-4256-bc2d-01f0e3136094
# ╟─8e1f15f2-8273-400d-97a1-86cd7d138a64
# ╟─e1d0a740-2c6b-4b72-ab3e-0c0fde78ed68
# ╟─b1b78dd1-112f-42b9-9d40-66123d64feac
# ╠═5d9e665a-7dd6-42cb-99f8-47ab69256feb
# ╟─3013da2b-e879-4730-a46f-ab3a94e6cbb8
# ╠═f4cf76b0-e5f6-4920-93a7-0be1f1041987
# ╟─57322645-33db-4c87-a5c8-51d69e8e74e9
# ╟─4843fc90-c81f-4810-8da8-8b1ce0694a55
# ╟─db97a178-677b-48f7-b15d-fbe64d374128
# ╟─2a8a27d1-55ff-4471-80c3-28a7c605b97e
# ╟─6cbcc482-d188-4732-8926-2cae312e2788
# ╠═394b491d-36b0-48a8-92d0-b2abb1ffdc01
# ╟─9e6b9683-458c-4796-94b0-590b7325c11b
# ╟─1ac0a8f6-6803-43bc-bc1e-8ac8d8b3ef1a
# ╟─56f09d0a-67c5-4109-bde0-b66296986a1c
# ╟─4f107e76-56a3-4b87-a76b-ec2179b082df
# ╟─47337c99-d8e1-4488-9959-77963d555e7d
# ╠═f613474b-8e87-4c3f-b352-570543e3990c
# ╟─53af53ab-a701-47cc-8979-14d6bef5e074
# ╟─4a6c39bf-04b0-4758-ba67-8648ed21817a
# ╟─d72010fc-cd62-4313-a3b2-aea4f5b277c2
# ╟─17d84f46-a294-4c12-b706-72ebc15a98e0
# ╠═d595178d-1edf-45d4-b97f-d8c98a8c88a4
# ╟─4da69a7c-1ece-485d-a353-cb944c900b26
# ╠═ba7cc1b5-e6bd-4b8b-adac-fa1636927d5e
# ╟─9a9975ef-2ec1-41af-bd6f-0a76756c100b
# ╟─e291e3fb-3301-44e9-8a54-842ccba1cf20
# ╟─263b843b-e75c-45a4-98a9-4f9111858ae6
# ╟─fba0bff1-e7a9-4793-94bd-968cd82fd926
# ╟─048fc935-3507-496b-b00d-09078784fe06
# ╟─1a8ed4b8-74ab-4775-b01b-24b637d9288c
# ╟─40940718-bf05-4f6e-8244-7716c77ca10c
# ╟─183a420c-cac8-47d1-b9ae-dff435f10858
# ╟─46991e0d-8a9c-487c-a157-72690cf2b489
# ╟─89819051-5d7b-47f7-ad1c-94362c4679fe
# ╟─a1de0642-7ef4-4d4b-ad3b-c2829cda3393
# ╟─bef48346-8bd5-4f24-9a1e-cc91f52318d1
# ╟─a82fb7fd-861d-4843-8cc6-18a33d9a41b4
# ╟─01efad8f-6aea-45ee-b084-9b68b35fb61b
# ╠═1beb011b-26bf-4718-85b1-c16747e184cd
# ╟─2b4ca2be-25dd-4f9c-a2c2-c4845c0f048a
# ╟─56f568bb-f68e-4c92-b1fa-89ed1a422c8a
# ╟─aa8a8293-eeb8-4327-bfe7-c3d131a506fe
# ╟─bbf4261d-7747-47e4-b69d-1d293c6172c0
# ╟─076b87bb-261d-4dca-a520-67f760a22d84
# ╠═b7bd9f6d-8df0-4cd4-b51c-2585da332d4c
# ╟─5a1b36d2-3065-4fb2-85e3-c5f131f2cd26
# ╟─4cea17b5-b132-4de2-9294-db3ebeccbbd3
# ╟─7cb3af49-e0c4-4692-9efd-605d53189a54
# ╟─34328fe5-caf3-441c-8957-0190ecc62170
# ╟─ef6598c7-2737-424a-8bbb-51298c9421a0
# ╟─17f6efeb-ccbb-4aa0-bf02-056195f91604
# ╟─9ee3d314-acee-4646-8139-6fd1706e7292
# ╟─de24c4e8-65c8-4952-a87f-5196ef74891a
# ╟─3a32d50b-75c5-407d-b0c4-2c4866ce7797
# ╟─db507c85-6847-4e29-aaf5-7314b22abe2d
# ╟─083fd2fa-582a-414f-87b3-012a3509a049
# ╟─71ea919d-fcf1-4282-8d71-2888406d31d0
# ╟─b3a6b457-b7ff-481f-9192-f421ba34dbf5
# ╟─bdce2e27-5fe4-47cc-b27e-7cfc2e320a45
# ╠═d21788c6-14d2-4237-b09e-0d869819aaf1
# ╟─249f483c-30ba-4ed3-9497-f774646999fc
# ╟─4ec999bf-a88f-40c8-92d6-0d578084ab5d
# ╟─67977543-3768-478f-b09f-a176a80f3f24
# ╟─f29b912a-7080-47f4-93cf-239a692bf337
# ╟─4e890dc7-db6a-4213-ac44-a2bb80e542c6
# ╟─ceaf76d6-d5ee-4334-a3d8-88373f5f0f31
# ╟─fb1e4d9c-cae8-4a42-b218-00f157ba7b60
# ╟─85794fff-8d0d-4ca3-bf94-b2aead8c9dd3
# ╟─7ad74395-20ad-4d02-92a9-788790382fa4
# ╠═6449b443-e4e4-40d9-aefc-a98d0dd65cea
# ╟─a0c5f1c5-4230-4630-8b4a-5515e7c49ebe
# ╟─1b72137a-1d39-49bc-9d2f-8d60262dccb1
# ╟─4a35b5c1-086f-4103-b8db-98886ddbbf9d
# ╟─2a6ad2c4-16f2-4a3c-9bf8-b198fe560626
# ╟─16ff57d8-3c79-4442-80f0-f8a43a949649
# ╟─04a47b4d-9a14-4225-a907-a63ac6ee70d6
# ╠═b84719e4-902e-4c13-a146-74c19c28ca5e
# ╟─e61790b0-ca31-44f8-bc2b-2bf610966400
# ╟─9c266479-0995-4f72-9942-0941bb0938fd
# ╟─1557cf00-e2db-49df-b1b7-07b83ff6e214
# ╟─0edaff9d-8630-435d-b2c5-4afea2ad7e5f
# ╟─3134e281-e748-4637-b182-97d27a14955f
# ╟─d09b1a7b-90fd-46aa-88b9-954dd01d8043
# ╟─1de083fd-d518-4c98-b27a-54ff74b6eb25
# ╟─a387c180-f6ef-4399-8df3-35c72e2330e2
# ╟─b382223f-9c78-4671-97ec-d898654dfd0a
# ╟─c03160d4-b063-4899-8df6-6ebd66cc4617
# ╠═3932efd6-8805-4b31-a724-0c45685fdb6f
# ╠═315db979-0472-45eb-a20a-9964821809f3
# ╠═279aea14-a829-4528-a7f6-d84eda6891ed
# ╠═87a8de33-2a99-4340-ad36-3d39abcc4de3
# ╠═14f1e18c-c521-4084-a13e-74a5e54fe572
# ╠═ba55ba2a-cb4e-4b31-b0df-9cf6593d56bc
# ╠═abf0c10e-b7d2-4a19-bfd5-2625015bccef
# ╟─ddb1b45e-df33-445c-b4b9-7611bc6f88a2
# ╟─658ced39-70ff-4fb1-8fdd-b314aca309ce
# ╟─8c4bf470-1408-4057-8d7e-216b1b38c1fa
# ╟─9222bbd2-a65e-425c-ac45-2abad5f43a5a
# ╟─a5d9fa3d-7715-4bc3-88ee-44871493c84d
# ╟─6c86adaf-2850-46f1-b421-e4ebc97e6631
# ╠═08b43888-cc2b-4b0c-88fd-7af121088db1
# ╠═4eb18bb0-5b04-11ef-0c2c-8747a3f06685
# ╟─ed7ac1ae-3da3-4a46-a34b-4b445d52a95f
# ╟─7b9ffd7c-3b93-4cfd-bed5-1590368ce987
# ╟─d0060e13-0aa0-495e-828e-084df48ef7e7
# ╟─d779340e-4dab-45c1-b8df-c0bcbae32a90
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
