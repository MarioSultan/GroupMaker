# Matrix Groups

In **GroupsMath**, matrix groups are represented by the `MatrixGroup` class. Unlike a `CayleyGroup`, which is defined by a complete Cayley table, a matrix group is defined directly by its matrix dimension, underlying field, and an optional membership condition.

The group operation is matrix multiplication.

---

## The `MatrixGroup` class

A matrix group is created using:

```python
MatrixGroup(n, field="R", condition=None, p=None)
```

where:

* `n` is the dimension of the square matrices.
* `field` specifies the underlying field.
* `condition` is an optional function imposing an additional condition.
* `p` is the prime defining the finite field when `field="F_p"`.

For example:

```python
from groupsmath import MatrixGroup

G = MatrixGroup(3, field="R")
```

represents the group of all invertible $3\times3$ real matrices, namely $GL_3(\mathbb{R})$.

By default, membership requires only that a matrix be invertible and have entries in the specified field.

---

## Supported fields

The `field` parameter currently accepts:

| Value | Mathematical set |
|---|---|
| `"Z"` | Integer matrices, over $\mathbb{Z}$ |
| `"R"` | Real matrices, over $\mathbb{R}$ |
| `"C"` | Complex matrices, over $\mathbb{C}$ |
| `"F_p"` | Matrices over the finite field $\mathbb{F}_p$ |

For example:

```python
MatrixGroup(2, field="Z")
MatrixGroup(2, field="R")
MatrixGroup(2, field="C")
MatrixGroup(2, field="F_p", p=5)
```

For convenience, a prime integer can also be passed directly: `G = MatrixGroup(2, 5)`, which is interpreted as a group over $\mathbb{F}_5$. The parameter `p` must be prime.

---

## Membership

The Python `in` operator checks whether a matrix belongs to a matrix group:

```python
A in G
```

GroupsMath verifies:

1. The matrix has shape $n\times n$.
2. Its entries belong to the selected field.
3. Its determinant is non-zero.
4. Any additional `condition` is satisfied.

For example:

```python
import numpy as np
from groupsmath import GL

G = GL(2)

A = np.array([[1, 2],
              [3, 4]])

A in G
# True
```

A singular matrix does not belong to `GL(2)`:

```python
A = np.array([[1, 2],
              [2, 4]])

A in G
# False
```

For $\mathbb{F}_p$, entries must be integer representatives in $\{0,1,\ldots,p-1\}$.

---

## Matrix operations

The group operation is available through `G.operation(A, B)`. For ordinary fields, GroupsMath computes $AB$.


The identity is obtained with `G.identity()`, and is the usual identity matrix $I_n$.

The inverse of a matrix is obtained with `G.inverse(A)`. For ordinary fields, GroupsMath uses the numerical matrix inverse. For finite fields, the inverse is computed modulo $p$ using:

$$
A^{-1}=(\det A)^{-1}\text{adj}(A)\pmod p.
$$

---

## Finite and infinite matrix groups

A matrix group is finite only when it is defined over a finite field:

```python
G = GL(2, 5)
```

Groups over `"Z"`, `"R"`, and `"C"` are treated as infinite.

The order is obtained using `G.order()`. Finite groups are enumerated and counted, while infinite groups return `float("inf")`

---

## Enumerating elements

Finite matrix groups can be enumerated with:

```python
G.elements()
```

This method is a generator.

For example:

```python
G = GL(2, 2)

for A in G.elements():
    print(A)
```

GroupsMath generates all possible $n\times n$ matrices over $\mathbb{F}_p$ and yields those satisfying the membership conditions.

Matrices are returned as tuples of tuples, for example:

```python
((1, 0),
 (0, 1))
```

Infinite matrix groups cannot be enumerated.

---

## Cayley tables

A Cayley table can only be generated for a finite matrix group:

```python
G.cayley_table()
```

Internally, GroupsMath converts the matrix group to a `CayleyGroup` and uses the usual Cayley table functionality.

Since the number of matrices grows rapidly with $n$ and $p$, this is mainly practical for small finite groups.

---

## Converting to a `CayleyGroup`

A finite matrix group can be converted using:

```python
C = G.toCayleyGroup()
```

The resulting object is a `CayleyGroup`, so its usual methods can be used:

```python
C.subgroups()
C.automorphisms()
C.automorphism_group()
```

Infinite matrix groups cannot be converted to `CayleyGroup`.

---

## Converting to an `ExplicitGroup`

Finite matrix groups can also be converted using:

```python
E = G.toExplicitGroup()
```

The resulting `ExplicitGroup` contains the matrix elements explicitly and uses matrix multiplication as its operation.

---

# Standard matrix groups

GroupsMath includes constructors for several important families of matrix groups.

---

## General linear groups: `GL`

```python
GL(n, field="R", p=None)
```

The general linear group is:

$$
GL_n(F)=\{A\in M_n(F)\mid\det(A)\neq0\}.
$$

For example:

```python
from groupsmath import GL

G = GL(3)
```

constructs $GL_3(\mathbb{R})$.

Over a finite field:

```python
G = GL(2, 5)
```

constructs $GL_2(\mathbb{F}_5)$.

---

## Special linear groups: `SL`

```python
SL(n, field="R", p=None)
```

The special linear group satisfies:

$$
SL_n(F)=\{A\in GL_n(F)\mid\det(A)=1\}.
$$

For example:

```python
from groupsmath import SL

G = SL(3)
G = SL(2, 5)
```

For finite fields, the determinant condition is interpreted modulo $p$.

---

## Orthogonal groups: `O`

```python
O(n, field="R")
```

Matrices in the orthogonal group satisfy:

$$
A^TA=I.
$$

For example:

```python
from groupsmath import O

G = O(3)
```

The current implementation does not support orthogonal groups over finite fields.

It also rejects even dimensions because a Witt index is not yet part of the constructor.

---

## Special orthogonal groups: `SO`

```python
SO(n, field="R")
```

Matrices satisfy:

$$
A^TA=I,
\qquad
\det(A)=1.
$$

For example:

```python
from groupsmath import SO

G = SO(3)
```

As with `O`, finite fields and even dimensions are currently not supported by this constructor.

---

## Unitary groups: `U`

```python
U(n)
```

Matrices satisfy:

$$
A^\dagger A=I,
$$

where $A^\dagger$ is the conjugate transpose.

For example:

```python
from groupsmath import U

G = U(3)
```

constructs $U(3)$.

---

## Special unitary groups: `SU`

```python
SU(n)
```

Matrices satisfy:

$$
A^\dagger A=I,
\qquad
\det(A)=1.
$$

For example:

```python
from groupsmath import SU

G = SU(2)
```

constructs $SU(2)$.

---

## Symplectic groups: `Sp`

```python
Sp(n)
```

The matrix dimension must be even.

The defining condition is:

$$
A^TJA=J,
$$

where:

$$
J=
\begin{pmatrix}
0&I\\
-I&0
\end{pmatrix}.
$$

For example:

```python
from groupsmath import Sp

G = Sp(4)
```

constructs the real symplectic group of $4\times4$ matrices, usually denoted $Sp_4(\mathbb{R})$.

An odd dimension raises a `ValueError`.

---

# The `MatrixSubgroup` class

Matrix subgroups are represented by:

```python
MatrixSubgroup(subgroup, group)
```

Both arguments must be instances of `MatrixGroup`.

The subgroup and its parent group must have:

* The same matrix dimension.
* The same field.
* The same value of `p` when working over $\mathbb{F}_p$.

A `MatrixSubgroup` inherits from both `MatrixGroup` and `Subgroup`.

---

## Custom matrix groups

The `condition` parameter can be used to define custom matrix groups.

The condition receives a NumPy matrix and must return `True` or `False`.

For example, invertible diagonal $3\times3$ real matrices can be represented by:

```python
import numpy as np
from groupsmath import MatrixGroup

G = MatrixGroup(
    3,
    field="R",
    condition=lambda A: np.allclose(A, np.diag(np.diag(A)))
)
```

A matrix belongs to this group only if it is invertible and satisfies the additional condition.

---

## Complete example

The following example constructs a finite matrix group:

```python
from groupsmath import GL

G = GL(2, 2)

print("Finite:", G.is_finite())
print("Order:", G.order())

for A in G.elements():
    print(A)

elements = list(G.elements())

A = elements[0]
B = elements[1]

print(G.operation(A, B))

C = G.toCayleyGroup()

print(C.order())
print(C.is_abelian())
```

The `MatrixGroup` class makes it possible to represent groups directly through matrices rather than through a precomputed Cayley table. Finite matrix groups can be enumerated and transformed into other group representations, while infinite matrix groups can still be studied through membership, multiplication, identity, inversion, and their defining conditions.
