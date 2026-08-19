# Elements

In **GroupMaker**, individual elements of a finite group can be represented explicitly using the `Element` class. While a `Group` object describes the complete algebraic structure through its Cayley table, an `Element` object represents one particular element together with the group to which it belongs.

The `Element` class allows group elements to be multiplied, raised to integer powers, inverted, and compared using a natural Python syntax.

---

## The `Element` class

An element is created using:

```python
Element(element, group)
```

where:

* `element` is the name of an element belonging to the group.
* `group` is an instance of the `Group` class.

For example:

```python
from groupmaker import Element
from groupmaker.precooked import C4

r = Element("$r$", C4)
```

The constructor checks that the provided element belongs to the group. If it does not, GroupMaker raises a `ValueError`. For example, `x = Element("x", C4)` raises an error because `"x"` is not an element of `C4`.

---

## Internal representation

When an `Element` object is created, GroupMaker stores both the visible element name and its internal index.

For example:

```python
r = Element("$r$", C4)

r.element
# '$r$'

r.index
# 1

r.group
# C4
```

The index is obtained from the position of the element inside `group.names`.

This distinction is important because the group operation is performed using the Cayley table, while the resulting `Element` object is created using the corresponding visible name.

---

## Multiplying elements

Two `Element` objects can be multiplied using the standard `*` operator:

```python
a * b
```

The multiplication is performed using the Cayley table of the group associated with the elements, therefore, both elements must belong to the same group. For example:

```python
r = Element("$r$", C4)

r * r
```

returns an `Element` representing: $r^2$.


The result is itself an `Element` object, so products can be chained:

```python
r * r * r
```

---

## Powers of elements

Elements can be raised to integer powers using the `**` operator:

```python
element ** n
```

For example:

```python
r = Element("$r$", C4)

r ** 2
```

represents: $r^2$.

Positive powers are calculated by repeated multiplication. The expression: `r ** 0` returns the identity element of the group. Negative powers are also supported and are constructed by repeatedly multiplying inverse elements $r^{-1}$.

---

## Comparing elements

Two `Element` objects can be compared using:

```python
a == b
```

Two elements are considered equal when both of the following are equal:

* Their element names.
* Their associated groups.

For example:

```python
r1 = Element("$r$", C4)
r2 = Element("$r$", C4)

r1 == r2
# True
```

The comparison does not only depend on the visible element name. The group is also part of the identity of an `Element` object.


---

## Complete example

The following example illustrates the main functionality of the `Element` class:

```python
from groupmaker import Element
from groupmaker.precooked import C4

# Create an element
r = Element("$r$", C4)

# Basic information
print(r)
print(r.element)
print(r.index)

# Products
print(r * r)
print(r * r * r)

# Powers
print(r ** 1)
print(r ** 2)
print(r ** 3)
print(r ** 0)
print(r ** -1)

# Inverse
rinv = r ** -1
print(rinv)

# Verification
print(r * rinv == Element("$e$", C4))
```

The `Element` class provides a direct object-oriented way of working with individual elements of a finite group. Instead of manually accessing entries of a Cayley table, group operations can be written using familiar algebraic notation such as: `a * b` and `a ** n` while the underlying computations are still performed entirely from the Cayley table of the associated `Group`.
