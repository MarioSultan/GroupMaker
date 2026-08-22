# Automorphisms

In **GroupsMath**, automorphisms are represented as permutations of the elements of a group. The library provides tools to test whether a permutation is an automorphism, compute all automorphisms of a finite group, construct the automorphism group $\text{Aut}(G)$, and represent automorphisms and actions between groups using the `Automorphism` and `AutomorphismFunction` classes.

---

## Automorphisms in GroupsMath

Mathematically, an automorphism of a group $G$ is an isomorphism from the group to itself:

$$
\varphi:G\longrightarrow G.
$$

Therefore, $\varphi$ must be bijective and preserve the group operation:

$$
\varphi(ab)=\varphi(a)\varphi(b)
$$

for every $a,b\in G$.

Since every finite group in GroupsMath is internally represented by indices

$$
0,1,\dots,n-1,
$$

an automorphism is represented by a tuple containing the image of every element.

For example, the tuple

```python
(0, 2, 1, 3)
```

represents the transformation

$$
0\mapsto0,\qquad
1\mapsto2,\qquad
2\mapsto1,\qquad
3\mapsto3.
$$

The tuple is interpreted as:

```python
phi[i] = image of i
```

The element names stored in `G.elements` are only used for displaying the group. Automorphisms themselves act on the internal indices of the Cayley table.

---

## Checking whether a permutation is an automorphism

The method

```python
G.is_automorphism(phi)
```

checks whether a tuple or list `phi` defines an automorphism of the group `G`.

The transformation must satisfy the following conditions:

1. It must be a permutation of all group indices.
2. It must preserve the identity element.
3. It must preserve the group operation.

In other words, GroupsMath verifies that

$$
\varphi(a\ast b)=\varphi(a)\ast\varphi(b)
$$

for every pair of elements of the group.

For example:

```python
from groupsmath.precooked import C4

phi = (0, 3, 2, 1)

C4.is_automorphism(phi)
# True
```

Here, `phi` represents the automorphism

$$
r\mapsto r^3,
$$

of the cyclic group $C_4$.

On the other hand:

```python
phi = (0, 1, 3, 2)

C4.is_automorphism(phi)
# False
```

because this permutation does not preserve the group operation.

If the argument is not a tuple or list, has the wrong length, or is not a complete permutation of the group elements, the method simply returns:

```python
False
```

---

## The `Automorphism` class

A valid automorphism can be explicitly represented by an instance of the `Automorphism` class.

The constructor is:

```python
Automorphism(phi, group)
```

where:

* `phi` is a tuple representing the permutation.
* `group` is the group on which the automorphism acts.

For example:

```python
from groupsmath import Automorphism
from groupsmath.precooked import C4

phi = (0, 3, 2, 1)

aut = Automorphism(phi, C4)
```

The constructor automatically checks whether `phi` is an automorphism of `C4`.

Therefore, the following will raise an error if the permutation does not define a valid automorphism:

```python
Automorphism((0, 1, 3, 2), C4)
# ValueError: The tuple phi has to be an automorphism of G.
```

The permutation defining an `Automorphism` object is stored in its `phi` attribute and the group associated with the automorphism is stored in its `group`attribute.


---

## Printing an automorphism

Printing an `Automorphism` object displays its permutation:

```python
print(aut)
# (0, 3, 2, 1)
```

Thus, the tuple representation can be directly inspected without accessing the `phi` attribute explicitly.

---

## Computing all automorphisms of a group

All automorphisms of a group can be computed using:

```python
G.automorphisms()
```

This method returns a list of `Automorphism` objects.

For example:

```python
from groupsmath.precooked import C4

auts = C4.automorphisms()

for aut in auts:
    print(aut)
```

The result is a list containing the automorphisms of $C_4$.

Each element of the returned list is an instance of `Automorphism`, so its permutation and associated group can be accessed normally:

```python
for aut in C4.automorphisms():
    print(aut.phi)
```

The number of automorphisms can therefore be obtained with:

```python
len(G.automorphisms())
```

Mathematically, this is the order of the automorphism group:

$$
|\text{Aut}(G)|.
$$


---

## The automorphism group: `G.automorphism_group()`

The set of all automorphisms of a group forms another group under composition:

$$
\text{Aut}(G)
=
\{\varphi:G\to G\mid\varphi\text{ is an automorphism}\}.
$$

In GroupsMath, this group can be constructed using:

```python
G.automorphism_group()
```

For example:

```python
from groupsmath.precooked import C4

AutC4 = C4.automorphism_group()
```

The result is a standard `CayleyGroup` object.

Therefore, all ordinary `CayleyGroup` methods can also be used with an automorphism group:

```python
AutC4.order()
AutC4.is_abelian()
AutC4.is_cyclic()
AutC4.cayley_table()
AutC4.subgroups()
```

This makes $\text{Aut}(G)$ behave computationally exactly like any other group created in GroupsMath.


---

## The `AutomorphismFunction` class

The `AutomorphismFunction` class represents a list of automorphisms associated with a group.

$$
f: \{0,1,...,n-1\}\rightarrow \text{Aut}(G)
$$

Its constructor is:

```python
AutomorphismFunction(function, group)
```

where:

* `function` is a list containing `Automorphism` objects or tuples representing automorphisms. The element `function[i]` represents $f(i)$, where the result is an automorphism.
* `group` is the group associated with the automorphisms.

For example:

```python
from groupsmath import AutomorphismFunction
from groupsmath.precooked import C4

f = AutomorphismFunction(
    [
        (0, 1, 2, 3),
        (0, 3, 2, 1),
        (0, 1, 2, 3),
        (0, 3, 2, 1)
    ],
    C4
)
```

The tuples are validated and stored internally as automorphisms.

The resulting list is available through:

```python
f.phi
```

which gives:

```python
[
    (0, 1, 2, 3),
    (0, 3, 2, 1),
    (0, 1, 2, 3),
    (0, 3, 2, 1)
]
```

An `AutomorphismFunction` object can also be created from `Automorphism` objects:

```python
id = Automorphism((0, 1, 2, 3), C4)
inv = Automorphism((0, 3, 2, 1), C4)

f = AutomorphismFunction([id, inv, id, inv], C4)
```

Both representations can be mixed inside the list.


---

## Important notes

* Automorphisms in GroupsMath are represented internally using indices, not the visible names in `G.elements`.

* The tuple `(0, 1, ..., n-1)` represents the identity automorphism.

* `G.automorphisms()` returns a list of `Automorphism` objects.

* `G.automorphism_group()` returns a standard `CayleyGroup` object representing $\text{Aut}(G)$.

* The group operation in $\text{Aut}(G)$ is composition of automorphisms.

* `AutomorphismFunction` stores a list of automorphisms and is primarily intended for constructions involving actions, such as semidirect products.

* When constructing a semidirect product, the position of each automorphism in `f.phi` is associated with the corresponding internal index of the acting group.

---

## Example: studying the automorphisms of a group

The following example computes all automorphisms of a group, constructs its automorphism group, and studies some of its properties:

```python
from groupsmath.precooked import C4

# Compute all automorphisms
auts = C4.automorphisms()

print("Automorphisms:")
for aut in auts:
    print(aut)

# Construct Aut(C4)
AutC4 = C4.automorphism_group()

print("Order:", AutC4.order())
print("Is cyclic?", AutC4.is_cyclic())
print("Is abelian?", AutC4.is_abelian())

# Display its Cayley table
AutC4.cayley_table()
```

Since `AutC4` is itself a `CayleyGroup`, it can be studied using the same methods as every other group in the library.

This illustrates a general idea used throughout GroupsMath: algebraic objects constructed from other algebraic objects can themselves be represented and manipulated using the same general `CayleyGroup` interface.
