# Cayley table visualization

GroupMaker allows you to visualize the Cayley table of any `Group` object using `matplotlib.pyplot`.

The Cayley table contains the result of the group operation for every ordered pair of elements. GroupMaker provides a dedicated method for displaying this information graphically.

---

## Displaying a Cayley table

The Cayley table of a group can be displayed using:

```python
G.cayley_table()
```

For example:

```python
from groupmaker.precooked import S3

S3.cayley_table()
```

The result is a graphical representation of the Cayley table, with the elements of the group shown along both axes.

By default, if the group has at most 20 elements, the result of each operation is also written inside its corresponding cell.

---

## Adding a title

A custom title can be added to the Cayley table using the `title` argument:

```python
S3.cayley_table(title="Cayley table of S3")
```

For example:

```python
C4.cayley_table(title="The cyclic group C4")
```

The title is displayed above the generated table.

---

## Changing the colormap

The appearance of the Cayley table can be modified using the `colormap` argument. Since the visualization is generated with `matplotlib`, any valid matplotlib colormap can be used. For example:

```python
S3.cayley_table(colormap="plasma")
```

or:

```python
S3.cayley_table(colormap="coolwarm")
```

GroupMaker also provides predefined colormaps, including: `rainbow`, `tgl` and `white`. The colormap `rainbow` is the default used by the method. Therefore, the following call:

```python
S3.cayley_table()
```

is equivalent to:

```python
S3.cayley_table(colormap=rainbow)
```

provided that the remaining arguments keep their default values.

---

## Showing element names inside the table

The `names` argument controls whether the names of the elements are displayed inside each cell.

To explicitly display them:

```python
S3.cayley_table(names=True)
```

To hide them:

```python
S3.cayley_table(names=False)
```

For example, the following produces a more compact visualization:

```python
S4.cayley_table(names=False)
```

This can be useful when working with larger groups.


> By default, **GroupMaker** automatically decides whether the element names should be written inside the cells. If the order of the group is at most $20$, the names are displayed. If the group has more than $20$ elements, they are hidden. This behavior can always be overridden manually writing `names=True` or `names=False`.

---

## Custom element names

The visualization always uses the names stored in `G.names` for the labels of the rows and columns. For example, consider:

```python
G = Group(
    [[0, 1], [1, 0]],
    names=["e", "a"]
)

G.cayley_table()
```

The axes and the values displayed inside the table will use `"e"` and `"a"` instead of the internal indices `0` and `1`.

This makes it possible to construct groups using numerical Cayley tables while displaying the elements using conventional mathematical notation. For example, a cyclic group may use:

```python
C5.names
# ['$e$', '$r$', '$r^{2}$', '$r^{3}$', '$r^{4}$']
```

and these names will automatically be used when plotting `C5.cayley_table()`.

---

## Complete examples

A simple Cayley table:

```python
from groupmaker.precooked import C4

C4.cayley_table()
```

A table with a custom title:

```python
from groupmaker.precooked import D4

D4.cayley_table(title="Symmetries of a square")
```

A table using a different colormap:

```python
from groupmaker.precooked import Q8

Q8.cayley_table(
    title="Quaternion group Q8",
    colormap="plasma"
)
```

A compact table without displaying the results inside each cell:

```python
from groupmaker.precooked import S4

S4.cayley_table(
    title="Cayley table of S4",
    names=False
)
```

The visualization method can therefore be adapted to both small groups, where displaying every element is useful, and larger groups, where a simpler graphical representation may be easier to read.
