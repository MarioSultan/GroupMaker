# Subgroups

Subgroups are one of the central structures in finite group theory. In GroupMaker, subgroups are represented by instances of the `Subgroup` class.

A `Subgroup` object inherits from the `Group` class, which means that all the usual methods available for groups can also be applied to subgroups.

This document explains how to obtain subgroups from a group and how to work with the additional functionality provided by the `Subgroup` class.


---

## Obtaining subgroups

All subgroups of a group $G$ can be obtained using:

```python
G.subgroups()
```

The returned list includes:

- The trivial subgroup $\{e\}$.
- All proper subgroups of $G$.
- The group $G$ itself.

For example:

```python
from groupmaker.precooked import *
S3.subgroups()
```

Each element of the returned list is an instance of the `Subgroup` class.

You can therefore iterate through them:

```python
for H in S3.subgroups():
    print(H)
```

Since a `Subgroup` behaves as a group, its usual properties can also be studied:

```python
for H in S3.subgroups():
    print(H.order())
    print(H.is_abelian())
```

---

## Proper subgroups

The proper subgroups of $G$ can be obtained using:

```python
G.proper_subgroups()
```

In GroupMaker, this method returns the non-trivial proper subgroups of the group.

For example:

```python
for H in S3.proper_subgroups():
    print(H)
```

The trivial subgroup and the whole group are not included in this result.

If all subgroups are required, including both of these cases, use `G.subgroups()`.

---

## The `Subgroup` class

A subgroup can also be created explicitly using the `Subgroup` constructor.

The constructor requires two `Group` objects:

```python
Subgroup(subgroup, group)
```

The first argument represents the subgroup, while the second argument represents the parent group.

For example:

```python
G = C4

H_group = Group(
    [[0, 1], [1, 0]],
    names=["$e$", "$r^{2}$"]
)

H = Subgroup(H_group, G)
```

The resulting object `H` is a subgroup of `G`.

When a `Subgroup` object is created, GroupMaker checks that:

- Both arguments are instances of `Group`.
- Every element of the proposed subgroup belongs to the parent group.
- The proposed subgroup actually forms a valid subgroup of the parent group.

If these conditions are not satisfied, an error is raised.

---

## Inherited group functionality

Since `Subgroup` inherits from `Group`, the usual group methods are available for every subgroup. For example:

```python
H.order()
H.identity()
H.element_orders()
H.order_distribution()
H.is_cyclic()
H.is_abelian()
H.center()
H.subgroups()
H.cayley_table()
```

This makes it possible to study a subgroup in exactly the same way as any other `Group` object. For example:

```python
H = C4.proper_subgroups()[0]

print(H.order())
print(H.is_cyclic())

H.cayley_table()
```

---

## The parent group

Every `Subgroup` object stores a reference to its parent group in the attribute:

```python
H.group
```

For example:

```python
H = S3.proper_subgroups()[0]

H.group     #S3
```

The subgroup itself is also stored internally as `H.subgroup`. The names of the elements of the parent group are available through `H.gnames` and the Cayley table of the parent group through `H.gcayley`.

These attributes are mainly used internally by GroupMaker to perform operations that depend on the relationship between a subgroup and its parent group.

---

## Cosets

Let $H$ be a subgroup of $G$ and let $g\in G$.

The left coset of $H$ determined by $g$ is

$$
gH=\{gh\mid h\in H\},
$$

while the right coset is

$$
Hg=\{hg\mid h\in H\}.
$$

In GroupMaker, cosets can be obtained using:

```python
H.coset(element)
```

By default, this returns the left coset and uses the element names of the parent group. For example:

```python
H = C4.proper_subgroups()[0]

H.coset("$r$")
```

The side can be specified explicitly with `H.coset("$r$", side="left")` or `H.coset("$r$", side="right")`.

The default behavior is:

```python
H.coset(element, side="left", return_names=True)
```

---

## Returning indices instead of names

By default, `coset()` returns the names of the elements. For example:

```python
H.coset("$r$")
```

To obtain the internal indices instead, use:

```python
H.coset("$r$", return_names=False)
```

This can be useful when the coset information is intended for further computation. For example:

```python
H.coset("$r$", side="left", return_names=False)
```

The argument `element` must correspond to one of the names stored in the parent group.

---

## Normal subgroups

A subgroup $H\leq G$ is normal if its left and right cosets coincide for every element of $G$:

$$
gH=Hg
$$

for every $g\in G$.

In GroupMaker, normality can be checked using:

```python
H.is_normal()
```

For example:

```python
for H in S3.subgroups():
    print(H.is_normal())
```

The normal subgroups of a group can also be obtained directly from the parent group:

```python
G.normal_subgroups()
```

For example:

```python
S3.normal_subgroups()
```

The result contains all subgroups of `G` for which `H.is_normal()` returns `True`.

---

## Quotient groups

A quotient group $G/H$ can only be constructed when $H$ is a normal subgroup of $G$.

If `H` is a `Subgroup` object, the quotient can be obtained using:

```python
G.quotient(H)
```

For example:

```python
Q = S3.quotient(H)
```

Finally, the `/` operator can be used:

```python
Q = S3 / H
```

Both approaches construct the quotient group associated with the normal subgroup `H`.

The resulting object is an instance of `Group`:

```python
Q.order()
Q.is_cyclic()
Q.is_abelian()
Q.cayley_table()
```

If `H` is not normal, GroupMaker raises an error when attempting to construct the quotient group.

---

## A complete example

The following example obtains the subgroups of $S_3$, identifies the normal ones, and constructs a quotient group.

```python
from groupmaker.precooked import S3

# Display all subgroups
for H in S3.subgroups():
    print(H)
    print("Order:", H.order())
    print("Normal:", H.is_normal())
    print()

# Select the normal subgroup of order 3
for H in S3.normal_subgroups():
    if H.order() == 3:
        A3 = H

# Construct S3 / A3
Q = S3 / A3

print(Q.order())        # 2
print(Q.is_cyclic())    # True

Q.cayley_table(title="Quotient group S3 / A3")
```

The `Subgroup` class therefore provides the connection between the internal structure of a subgroup and the larger group to which it belongs, while retaining all the functionality of the `Group` class.
