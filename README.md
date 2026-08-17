# GroupMaker – v0.1.0

![GroupMaker Logo](https://raw.githubusercontent.com/MarioSultan/GroupMaker/main/logo.png)


**GroupMaker** is a Python module for constructing, studying, and manipulating finite groups.

The project aims to provide a simple and intuitive interface for working with finite groups and their algebraic properties.

> **Note:** GroupMaker is currently under active development. Its API is still evolving and may change before version `1.0.0`.

---

## Installation

GroupMaker is currently distributed through GitHub and has not yet been published on PyPI.

To install the latest development version directly from GitHub, run:

```bash
pip install git+https://github.com/MarioSultan/GroupMaker.git
```

Alternatively, you can visit https://github.com/MarioSultan/GroupMaker and download all the files.


---

## Quick start

GroupMaker represents finite groups using the `Group` class. A group can then be constructed using the `Group` constructor:

```python
from groupmaker import *

G = Group(...)
```

Once a group has been created, its properties and operations can be accessed through its methods.

For example, the order of a group can be obtained with:

```python
G.order()
```

and its Cayley table with:

```python
G.cayley_table()
```

> The API is currently being migrated from a function-based interface to an object-oriented interface based on the `Group` class. Some functionality may therefore change during development.

---

## Features

GroupMaker is designed to provide tools for working with finite groups, including:

* Construction and representation of finite groups.
* Computation of Cayley tables.
* Determination of the order of a group.
* Study of subgroups.
* Detection of algebraic properties such as cyclicity and commutativity.
* Computation and exploration of other structural properties of finite groups.

The available functionality will expand as the project develops.

---

## Documentation

The documentation is currently being developed alongside the new object-oriented API.

Detailed documentation will be available in the [`docs`](docs/) directory.

---

## Project status

GroupMaker is currently under active development.

The project is undergoing a transition from a collection of functions to an object-oriented interface centred around the `Group` class.

The API should therefore be considered **unstable** until the `1.0.0` release.

Breaking changes may occur between development versions.

---

## Contributing

Contributions, suggestions, and bug reports are welcome.

If you find a bug or have an idea for improving GroupMaker, please open an issue in the GitHub repository.
