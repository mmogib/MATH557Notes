### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> title = "TETM241 | MATH557 | NOTES"
#> 
#>     [[frontmatter.author]]
#>     name = "Dr. Mohammed Alshahrani"
#>     url = "https://mshahrani.website/"

using Markdown
using InteractiveUtils

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

# ╔═╡ 85794fff-8d0d-4ca3-bf94-b2aead8c9dd3
TableOfContents(title="📚 MATH557: Applied Linear Algebra", indent=true,depth=4)

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

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "ed431b402d73b8c6440b3ee125486d75a4458251"

[[deps.ADTypes]]
git-tree-sha1 = "99a6f5d0ce1c7c6afdb759daa30226f71c54f6b0"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.7.1"

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
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "f61b15be1d76846c0ce31d3fcfac5380ae53db6a"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.37"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
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
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "f54c23a5d304fb87110de62bace7777d59088c34"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.15.0"

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

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "95f5c7e2d177b7ba1a240b0518038b975d72a8c0"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.7"

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
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
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
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

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

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "0e0a1264b0942f1f3abb2b30891f2a590cc652ac"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.110"

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

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

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
git-tree-sha1 = "94997910aca72897524d2237c41eb852153b0f65"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
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

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "3f74912a156096bd8fdbef211eff66ab446e7297"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+0"

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

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "801aef8228f7f04972e596b09d4dba481807c913"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.4"

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
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

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
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

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
git-tree-sha1 = "2787db24f4e03daf859c6509ff87764e4182f7d1"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.16"
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
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

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
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

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
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

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
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

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
git-tree-sha1 = "27d162f37cc29de047b527dab11a826dd3a650ad"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

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

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "5c1d1d9361e1417e5a065e1f84dc3686cbdaea21"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "d0a6b1096b584a2b88efb70a92f8cb8c881eb38a"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.6"

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
version = "0.3.23+4"

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
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

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
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

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
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

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
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

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

[[deps.PlutoDevMacros]]
deps = ["AbstractPlutoDingetjes", "DocStringExtensions", "HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Pkg", "Random", "TOML"]
git-tree-sha1 = "c3839362a712e6d9c2845d179edafe74371cb77b"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.7.4"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoDevMacros", "PlutoUI", "REPL"]
git-tree-sha1 = "93d8c75734da9192d0639406fe6fb446be0fba4f"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.12"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "d7f3f63331c7c8c81245b4ee2815b7d496365833"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.23"

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

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

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
git-tree-sha1 = "e237232771fdafbae3db5c31275303e056afaa9f"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.10.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e60724fd3beea548353984dc61c943ecddb0e29a"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.3+0"

[[deps.Roots]]
deps = ["Accessors", "ChainRulesCore", "CommonSolve", "Printf"]
git-tree-sha1 = "3484138c9fa4296a0cf46a74ca3f97b59d12b1d0"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.1.6"

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
git-tree-sha1 = "2803cab51702db743f3fda07dd1745aadfbf43bd"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.5.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "7f0e208db50f5fee2386b6d8dc9608d580059331"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.48.1"

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
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "23b02c588ac9a17ecb276cc62ab37f3e4fe37b32"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.9"
weakdeps = ["SparseArrays"]

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "20ad3e7c137156c50c93c888d0f2bc5b7883c729"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.4.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

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
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

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

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

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
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

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
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
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

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "c9fce29fb41a10677e24f74421ebe31220b81ad0"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.28"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fb099adbd7504f1e68b4512828e9d94197a8b889"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "9345b7b8a2923abaf6089d9f7306bb712fecb840"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "2.1.2"

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

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "5a13ae8a41237cff5ecf34f73eb1b8f42fff6531"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.24"

[[deps.TranscodingStreams]]
git-tree-sha1 = "e84b3a11b9bece70d14cce63406bbc79ed3464d2"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.2"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

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
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

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
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

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
version = "1.52.0+1"

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
# ╟─c59dffb4-c419-461b-8096-e27171be0a87
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
# ╠═6d988a2c-f6ae-418d-9cd4-11207ec7a2ce
# ╠═02c7ec2d-e57f-49d5-a8c3-a6f0637e026f
# ╟─27e009e2-2e77-4e58-a2d8-50a9efb99e32
# ╟─9579da39-7241-4266-918a-0c46b656cbf3
# ╟─6e8ecdb1-7d32-429f-9d23-8b79db0cc75d
# ╟─85794fff-8d0d-4ca3-bf94-b2aead8c9dd3
# ╠═4eb18bb0-5b04-11ef-0c2c-8747a3f06685
# ╟─ed7ac1ae-3da3-4a46-a34b-4b445d52a95f
# ╟─7b9ffd7c-3b93-4cfd-bed5-1590368ce987
# ╟─d0060e13-0aa0-495e-828e-084df48ef7e7
# ╟─d779340e-4dab-45c1-b8df-c0bcbae32a90
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
