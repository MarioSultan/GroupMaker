# Products of groups

In **GroupMaker**, new groups can be constructed from existing ones using group product constructions. The library supports **direct products**, **direct powers**, and **semidirect products**.

Because product groups in GroupMaker return standard instances of the `Group` class, all regular group properties and methods—such as `.order()`, `.is_abelian()`, `.subgroups()`, and `.cayley_table()`—apply to them directly.

---

## Direct product of two groups

Given two finite groups $A$ and $B$, their direct product $A \times B$ is the Cartesian product set $A \times B$ with the componentwise operation:

$$
(a_1, b_1) \ast (a_2, b_2) = (a_1 a_2, b_1 b_2)
$$

In GroupMaker, the direct product can be constructed using the function `direct_product(A, B)` or using the standard multiplication operator `*`:

```python
from groupmaker import direct_product
from groupmaker.precooked import C2, C3

# Using the function
G = direct_product(C2, C3)

# Using the * operator
G = C2 * C3
```

The resulting group has order $|A| \times |B|$.

### Element names in direct products

The names of the elements in $A \times B$ are formed by joining the corresponding element names of $A$ and $B$ with a comma:

```python
from groupmaker.precooked import C2, C3

G = C2 * C3
print(G.names)
# ['$e$,$e$', '$e$,$r$', '$r$,$e$', '$r$,$r$']
```

---

## Direct powers of a group

The direct power $G^n$ represents the direct product of $n$ copies of $G$:

$$
G^n = \underbrace{G \times G \times \dots \times G}_{n \text{ times}}
$$

In GroupMaker, direct powers can be constructed using the function `direct_power(G, n)` or using the exponentiation operator `**`:

```python
from groupmaker import direct_power
from groupmaker.precooked import C2

# Using the function
G3 = direct_power(C2, 3)

# Using the ** operator
G3 = C2 ** 3
```

The integer exponent $n$ must be at least $2$. Attempting to compute `G ** 1` or `G ** 0` will raise a `ValueError`.

For example, the Klein four-group $V_4$ is isomorphic to $C_2 \times C_2$:

```python
from groupmaker.precooked import C2

V4_alt = C2 ** 2

print("Order:", V4_alt.order())          # 4
print("Is abelian?", V4_alt.is_abelian()) # True
print("Is cyclic?", V4_alt.is_cyclic())   # False
```

---

## Semidirect products

Let $A$ and $B$ be two groups, and let $\phi: B \to \text{Aut}(A)$ be a group homomorphism from $B$ into the automorphism group of $A$. 

The **semidirect product** $A \rtimes_\phi B$ is the group whose underlying set is $A \times B$, with the group operation defined by:

$$
(a_1, b_1) \ast (a_2, b_2) = (a_1 \cdot \phi(b_1)(a_2), \, b_1 \cdot b_2)
$$

Here, $A$ is a normal subgroup of $A \rtimes_\phi B$.

### Construction in GroupMaker: `semidirect_product(A, B, f)`

To construct a semidirect product in GroupMaker, you must specify:

1. `A`: The normal group $A$.
2. `B`: The acting group $B$.
3. `f`: An `AutomorphismFunction` instance defining the mapping $\phi: B \to \text{Aut}(A)$.

```python
semidirect_product(A, B, f)
```

The object `f` maps each internal index $b_1 \in \{0, 1, \dots, |B|-1\}$ of group $B$ to an automorphism of group $A$, represented as a tuple or an `Automorphism` instance.

---

## Complete Example: Constructing the Dihedral Group $D_3$ as $C_3 \rtimes C_2$

The dihedral group $D_3$ (of order 6) can be constructed as a semidirect product $C_3 \rtimes_\phi C_2$, where $C_2 = \{0, 1\}$ acts on $C_3$ by inversion:

* For element $0 \in C_2$ (identity): map to the identity automorphism $(0, 1, 2)$ of $C_3$.
* For element $1 \in C_2$ (generator): map to the inversion automorphism $(0, 2, 1)$ of $C_3$.

```python
from groupmaker import semidirect_product, AutomorphismFunction
from groupmaker.precooked import C3, C2

# Define automorphisms of C3
id_aut = (0, 1, 2)   # Identity automorphism: r -> r
inv_aut = (0, 2, 1)  # Inversion automorphism: r -> r^2

# Map C2 elements [0, 1] to automorphisms of C3
f = AutomorphismFunction([id_aut, inv_aut], C3)

# Construct C3 semidirect C2
D3_constructed = semidirect_product(C3, C2, f)

# Inspect properties
print("Order:", D3_constructed.order())          # 6
print("Is abelian?", D3_constructed.is_abelian()) # False
print("Element names:", D3_constructed.names)
# ['$e$,$e$', '$e$,$r$', '$r$,$e$', '$r$,$r$', '$r^{2}$,$e$', '$r^{2}$,$r$']
```

---

## Properties of Product Groups

Since products in GroupMaker return standard `Group` instances, all analysis tools apply seamlessly:

```python
from groupmaker.precooked import C2, C3

G = C2 * C3

# Basic properties
print("Order:", G.order())                  # 6
print("Is cyclic?", G.is_cyclic())          # True
print("Order distribution:", G.order_distribution()) # {1: 1, 2: 1, 3: 2, 6: 2}

# Subgroups and Cayley table
print("Number of subgroups:", len(G.subgroups()))
G.cayley_table(title="Cayley Table of C2 x C3")
```
