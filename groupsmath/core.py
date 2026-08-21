# ———————————————————————————————————————————————————————————— #
#                                                              #
#                         GroupsMath                           #
#           A Python library for finite group theory           #
#                                                              #
#                   by: Mario Sultan Romero                    #
#                                                              #
# ———————————————————————————————————————————————————————————— #


#################### IMPORTS ####################

import  numpy               as      np
import  matplotlib.pyplot   as      plt
import  matplotlib.colors   as      mcolors
from    itertools           import  permutations, combinations, product
from    collections         import  Counter
from    math                import  gcd
from    abc                 import  abstractmethod, ABC


#################### DEFINITIONS ####################

__version__ = "0.4.1"

tgl_color = "#2A646E"
white = mcolors.LinearSegmentedColormap.from_list("white", ["white", "white"])
tgl = mcolors.LinearSegmentedColormap.from_list("tgl", ["white", tgl_color])
rainbow = mcolors.LinearSegmentedColormap.from_list("rainbow", ["#FF0000","#FF9100","#F2DA00","#2CDB00","#00DAE9","#1869FF"])#,"#7648FF","#D12BFF"])
rainbow8 = mcolors.LinearSegmentedColormap.from_list("rainbow8", ["#FF0000","#FF9100","#F2DA00","#2CDB00","#00DAE9","#1869FF","#7648FF","#FF1FF4"])

def info():
    print("———————————————————————————————————————————————————————————————————————————")
    print("                                                                           ")
    print("                          > GroupsMath v"+__version__+"                              ")
    print("              The educational Python library for Group Theory              ")
    print("                                                                           ")
    print("GroupsMath is a Python module for constructing, studying, and manipulating ")
    print("finite groups. The project aims to provide a simple and intuitive interface")
    print("for working with finite groups and their algebraic properties.             ")
    print("                                                                           ")
    print("You can view all the information of the proyect in PyPI and GitHub:        ")
    print("· PyPI https://pypi.org/project/GroupsMath/                                ")
    print("· GitHub https://github.com/MarioSultan/GroupsMath                         ")
    print("                                                                           ")
    print("The documentation of this library is in the following link:                ")
    print("· Documentation https://github.com/MarioSultan/GroupsMath/tree/main/docs   ")
    print("                                                                           ")
    print("———————————————————————————————————————————————————————————————————————————")


#################### GROUP CLASS ####################

class Group(ABC):
    @classmethod
    @abstractmethod
    def operation(self, a, b):
        pass

    @classmethod
    @abstractmethod
    def identity(self):
        pass

    @classmethod
    @abstractmethod
    def inverse(self, a):
        pass

class ExplicitGroup(Group):

    def __init__(self, elements, function, _skip_validation=False):

        if type(elements)!=list:
            raise ValueError("The first argument of a ExplicitGroup instance must be a list")

        if not _skip_validation:
            valid, message = _check_explicit_group(elements,function)
            if not valid:
                raise ValueError(message)

        self.elements = elements
        self.function = function

    #–> DUNDERS 

    def __str__(self):
        return str(self.toCayleyGroup())
    
    def __len__(self):
        return len(self.elements)
    
    def __contains__(self, element):
        return element in self.elements
    
    def __mul__(self, other):
        return direct_product(self.toCayleyGroup(),other.toCayleyGroup()).toExplicitGroup()
    
    def __pow__(self, n):
        return direct_power(self.toCayleyGroup(),n).toExplicitGroup()
    
    def __truediv__(self, subgroup):
        return self.toCayleyGroup().quotient(subgroup.toCayleyGroup()).toExplicitGroup()

    #–> TRANSFORMING METHODS 

    def toCayleyGroup(self):
        idx = {elem: i for i, elem in enumerate(self.elements)}
        cayley = [[idx[self.function(a, b)] for b in self.elements] for a in self.elements]
        return CayleyGroup(cayley, self.elements, _skip_validation=True)

    #–> GROUP METHODS 
    
    def operation(self, a, b):
        return self.function(self.elements[a],self.elements[b])

    def identity(self):
        neutro = None
        for e in range(self.order()):
            ok = True
            for a in range(self.order()):
                if self.function(e,a) != a:
                    ok = False
                    break
                if self.function(e,a) != a:
                    ok = False
                    break
            if ok:
                neutro = e
                break
        return neutro

    def inverse(self, a):
        e = self.identity()
        inv = None
        for i in self.order():
            if self.function(a,i)==e:
                inv = i
                break
        return inv

    #–> BASIC METHODS 

    def order(self):
        return len(self.elements)

    def element_orders(self):
        O = []
        for i in range(len(self.elements)):
            ik = self.function(self.elements[i],self.elements[i])
            o = 1
            while i!=ik:
                ik = self.function(self.elements[ik],self.elements[i])
                o+=1
            O.append(o)
        return O

    def order_distribution(self):
        l = sorted(self.element_orders())
        return dict(Counter(l))

    def is_cyclic(self):
        return self.order() in self.element_orders()

    def center(self):
        Z = []
        for i in self.elements:
            r = True
            for j in self.elements:
                if self.function(i,j)!=self.function(j,i):
                    r = False
                    break
            if r:
                Z.append(i)
        return(Z)

    def is_abelian(self):
        return self.order()==len(self.center())

    def cayley_table(self, title="", colormap=rainbow, names=None):
        return self.toCayleyGroup().cayley_table(self, title=title, colormap=colormap, names=names)

    #–> SUBGROUPS 

    def proper_subgroups(self):             # PROGRAMAR
        return self.toCayleyGroup().proper_subgroups()

    def subgroups(self):                    # PROGRAMAR
        return self.toCayleyGroup().subgroups()

    def normal_subgroups(self):             # PROGRAMAR
        return self.toCayleyGroup().normal_subgroups()

    def is_simple(self):                    # PROGRAMAR
        return self.toCayleyGroup().is_simple()

    def quotient(self, subgroup):           # PROGRAMAR
        return self.toCayleyGroup().quotient(subgroup.toCayleyGroup())

    #–> AUTOMORPHISMS 

    def is_automorphism(self, phi):
        return self.toCayleyGroup().is_automorphism(phi)

    def automorphisms(self):
        return self.toCayleyGroup().automorphisms()
    
    def automorphism_group(self):           # REVISAR (Aunque seguramente se quede así)
        return self.toCayleyGroup().automorphism_group().toExplicitGroup()

class CayleyGroup(Group):

    def __init__(self, cayley, elements=None, _skip_validation=False):

        if not _skip_validation:
            valid, message = _check_cayley_group(cayley)
            if not valid:
                raise ValueError(message)

        if elements is None:
            elements = list(range(len(cayley)))

        if len(elements) != len(cayley):
            raise ValueError("The number of elements must match the group order")

        self.cayley = cayley
        self.elements = elements

    #–> DUNDERS 

    def __str__(self):
        G = []
        for i in range(self.order()):
            g = []
            for j in range(self.order()):
                g.append(self.elements[self.cayley[i][j]])
            G.append(g)
        s = str(G)
        return f"{G}"

    def __len__(self):
        return self.order()

    def __contains__(self, element):
        return element in self.elements

    def __mul__(self, other):
        return direct_product(self,other)

    def __pow__(self, n):
        return direct_power(self,n)

    def __truediv__(self, subgroup):
        return self.quotient(subgroup)

    #–> TRANSFORMING METHODS 

    def toExplicitGroup(self):
        idx = {elem: i for i, elem in enumerate(self.elements)}
        def op(a, b):
            return self.elements[self.cayley[idx[a]][idx[b]]]
        return ExplicitGroup(list(self.elements), op, _skip_validation=True)
    
    #–> GROUP METHODS 

    def operation(self, a, b):
        return self.cayley[a][b]

    def identity(self):
            neutro = None
            for e in range(self.order()):
                ok = True
                for a in range(self.order()):
                    if self.cayley[e][a] != a:
                        ok = False
                        break
                    if self.cayley[a][e] != a:
                        ok = False
                        break
                if ok:
                    neutro = e
                    break
            return neutro

    def inverse(self, a):
        e = self.identity()
        inv = None
        for i in range(len(self.cayley)):
            if self.cayley[a][i]==e:
                inv = i
                break
        return inv

    #–> BASIC METHODS 

    def _print_group(self):
        print(f"CayleyGroup({self.cayley},{self.elements})")

    def order(self):
        return len(self.cayley)

    def element_orders(self):
        O = []
        for i in range(len(self.cayley)):
            ik = self.cayley[i][i]
            o = 1
            while i!=ik:
                ik = self.cayley[ik][i]
                o+=1
            O.append(o)
        return O

    def order_distribution(self):
        l = sorted(self.element_orders())
        return dict(Counter(l))

    def is_cyclic(self):
        return self.order() in self.element_orders()

    def center(self):
        G = self.cayley
        Z = []
        for i in range(len(G)):
            r = True
            for j in range(len(G)):
                if G[i][j]!=G[j][i]:
                    r = False
                    break
            if r:
                Z.append(i)
        return(Z)

    def is_abelian(self):
        return self.order()==len(self.center())

    def cayley_table(self, title="", colormap=rainbow, names=None):

        if names==None:
            if len(self.cayley)<=20:
                names=True
            else:
                names=False


        elements = range(len(self.cayley))
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(self.cayley, cmap=colormap)
        ax.xaxis.tick_top()

        names_in_axis = [f"${i}$" for i in self.elements]

        if names:
            ax.set_xticks(np.arange(len(elements)))
            ax.set_yticks(np.arange(len(elements)))
            ax.set_xticklabels(names_in_axis)
            ax.set_yticklabels(names_in_axis)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        if names:
            for i in range(len(elements)):
                for j in range(len(elements)):
                    color_texto = "black"
                    ax.text(
                        i,
                        j,
                        f"${self.elements[self.cayley[j][i]]}$",
                        ha="center",
                        va="center",
                        color=color_texto,
                        fontsize=12,
                    )
        
        plt.title(title)
        plt.tight_layout()
        plt.show()

    def delete_names(self):
        self.elements = [i for i in range(len(self.elements))]

    #–> SUBGROUPS 

    def proper_subgroups(self):
        G = self.cayley
        E = self.elements
        L = []
        
        for idxs in _graded_power_set_with_id(G):
            # 1. Comprobación rápida de clausura
            if _is_closed_subset(G, idxs):
                # 2. Construcción rápida omitiendo axiomas heredados
                sub_names = [E[j] for j in idxs]
                sub_matrix = _reset_renaming(_subset(G, idxs))
                
                sub_grp = CayleyGroup(sub_matrix, sub_names, _skip_validation=True)
                L.append(CayleySubgroup(sub_grp, self))
                
        return L

    def subgroups(self):
        # 1. Subgrupo trivial {e}
        e_name = self.elements[0]
        trivial_grp = CayleyGroup([[0]], [e_name])
        trivial_subgroup = CayleySubgroup(trivial_grp, self)

        # 2. Subgrupos propios
        subs = self.proper_subgroups()

        # 3. Subgrupo total (el propio grupo G)
        total_grp = CayleyGroup(self.cayley, self.elements)
        total_subgroup = CayleySubgroup(total_grp, self)

        return [trivial_subgroup] + subs + [total_subgroup]

    def normal_subgroups(self):
        return [sub for sub in self.subgroups() if sub.is_normal()]

    def is_simple(self):
        normals = self.normal_subgroups()
        return len(normals) == 2

    def quotient(self, subgroup):
        if not isinstance(subgroup, CayleySubgroup):
            raise TypeError("Argument must be an instance of CayleySubgroup")
        if subgroup.group is not self:
            raise ValueError("The subgroup does not belong to this group")
        return subgroup.quotient()

    #–> AUTOMORPHISMS 

    def is_automorphism(self, phi):

        if not isinstance(phi, (tuple, list)):
            return False

        n = self.order()
        if len(phi) != n:
            return False

        # 1. Biyectividad (debe ser una permutación completa de los índices 0..n-1)
        if set(phi) != set(range(n)):
            return False

        # 2. Conservación del elemento neutro (el índice 0)
        if phi[0] != 0:
            return False

        # 3. Preservación de la operación: phi(a * b) == phi(a) * phi(b)
        G = self.cayley
        for i in range(n):
            for j in range(n):
                if phi[G[i][j]] != G[phi[i]][phi[j]]:
                    return False

        return True

    def automorphisms(self):
        n = self.order()
        auts = []
        orders = self.element_orders()[:]
        for p in _order_preserving_permutations(orders):
            try:
                auts.append(Automorphism(p,self))
            except:
                pass

        return auts
    
    def automorphism_group(self):
        auts = [i.phi for i in self.automorphisms()]
        n = len(auts)
        dict_auts = {a: i for i, a in enumerate(auts)}

        # Tabla de Cayley de Aut(G) mediante composición de permutaciones: (a o b)[x] = a[b[x]]
        cayley_aut = []
        for i in range(n):
            row = []
            a = auts[i]
            for j in range(n):
                b = auts[j]
                comp = tuple(a[b[k]] for k in range(len(b)))
                row.append(dict_auts[comp])
            cayley_aut.append(row)

        names_aut = [str(auts[i]) for i in range(n)]
        return CayleyGroup(cayley_aut, names_aut, _skip_validation=True)

#class MatrixGroup(Group)

class Subgroup(ABC):
    pass

class CayleySubgroup(CayleyGroup, Subgroup):

    def __init__(self, subgroup:CayleyGroup, group:CayleyGroup):

        if not isinstance(subgroup, CayleyGroup) or not isinstance(group, CayleyGroup):
            raise ValueError("Both arguments must be instances of CayleyGroup")

        try:
            subgroup_indices = [group.elements.index(name) for name in subgroup.elements]
        except ValueError:
            raise ValueError("All elements of subgroup must belong to group")

        if not _form_subgroup(group.cayley, subgroup_indices):
            raise ValueError("The provided group is not a valid subgroup of the main group")

        # Todo lo que se puede hacer con grupos ahora también con subgrupos.
        super().__init__(subgroup.cayley, subgroup.elements)

        self.subgroup = subgroup
        self.group = group
        self.cayley = subgroup.cayley
        self.elements = subgroup.elements
        self.gcayley = group.cayley
        self.gelements = group.elements
        self._indices = subgroup_indices

    #–> DUNDERS 

    def __truediv__(self, other):
        """Permite usar la sintaxis natural G / H o H / H."""
        return self.quotient()

    def __le__(self, group):
        return group==self.group

    def __lt__(self, group):
        return group==self.group and self.subgroup.order()<self.group.order()

    #–> QUOTIENTS 

    def coset(self, element, side="left", return_names=True):
        if side not in ("left", "right"):
            raise ValueError("side must be either 'left' or 'right'")

        # Buscar directamente la primera aparición del elemento en los nombres del grupo padre
        try:
            elem_idx = self.gelements.index(element)
        except ValueError:
            raise ValueError(f"Element '{element}' is not present in parent group names.")

        coset_indices = []
        for h in self._indices:
            ah = self.gcayley[elem_idx][h] if side == "left" else self.gcayley[h][elem_idx]
            if ah not in coset_indices:
                coset_indices.append(ah)

        coset_indices.sort()

        if return_names:
            return [self.gelements[i] for i in coset_indices]
        return coset_indices

    def is_normal(self):
        for name in self.gelements:
            left = self.coset(name, side="left", return_names=False)
            right = self.coset(name, side="right", return_names=False)
            if left != right:
                return False
        return True

    def quotient(self):
        """Calcula el grupo cociente G/H devolviendo una instancia de CayleyGroup."""
        if not self.is_normal():
            raise ValueError("The subgroup must be normal to construct a quotient group.")

        # 1. Obtener todas las clases laterales (cosets) únicas expresadas como listas de nombres
        cosets = []
        for name in self.gelements:
            c = self.coset(name, side="left", return_names=True)
            if c not in cosets:
                cosets.append(c)

        n_cosets = len(cosets)

        # 2. Mapear cada elemento de G al índice de su clase lateral en 'cosets'
        elem_to_coset = {}
        for idx, c in enumerate(cosets):
            for name in c:
                elem_to_coset[name] = idx

        # 3. Construir la tabla de Cayley del grupo cociente
        quotient_cayley = []
        for i, c1 in enumerate(cosets):
            row = []
            rep1_idx = self.gelements.index(c1[0])
            for j, c2 in enumerate(cosets):
                rep2_idx = self.gelements.index(c2[0])
                # Producto en el grupo padre: g1 * g2
                prod_idx = self.gcayley[rep1_idx][rep2_idx]
                prod_name = self.gelements[prod_idx]
                row.append(elem_to_coset[prod_name])
            quotient_cayley.append(row)

        # 4. Asignar nombres representativos a los cosets, p. ej. "{e, r}" o "gH"
        quotient_names = [f"{{{','.join(str(e) for e in c)}}}" for c in cosets]

        return CayleyGroup(quotient_cayley, quotient_names, _skip_validation=True)

#class ExplicitSubgroup(ExplicitGroup, Subgroup)

#class MatrixSubgroup(MatrixGroup, Subgroup)

class Element:

    def __init__(self,element,group):

        if not isinstance(group, CayleyGroup):
            raise TypeError("Argument must be an instance of CayleyGroup")

        if not element in group:
            raise ValueError("The element must belong to the group.")

        self.element = element
        self.group = group
        self.index = group.elements.index(element)

    #–> DUNDERS 

    def __str__(self):
        return self.element

    def __mul__(self, other):
        return Element(self.group.elements[self.group.cayley[self.index][other.index]],self.group)

    def __pow__(self,k):
        if type(k)!=int:
            raise ValueError("n must be an integer")
        if k==1:
            return self
        elif k>1:
            p = self
            for i in range(k-1):
                p = p * self
            return p
        elif k==0:
            return Element(self.group.identity(),self.group)
        elif k==-1:
            return self.inverse()
        elif k<-1:
            p = self.inverse()
            for i in range(k-1):
                p = p * self
            return p

    def __eq__(self,other):
        return self.element==other.element and self.group==other.group

    #–> OTHER FUNCTIONS 

    def inverse(self):
        inv = None
        for i in self.group.elements:
            if self*Element(i,self.group)==Element(self.group.elements[self.group.identity()],self.group):
                inv = i
                break
        return Element(inv,self.group)

class Automorphism:

    def __init__(self,phi:tuple,group:CayleyGroup):

        if not group.is_automorphism(phi):
            raise ValueError("The tuple phi has to be an automorphism of G.")

        self.phi = phi
        self.group = group

    def __len__(self):
        return len(self.phi)

    def __str__(self):
        return str(self.phi)

class AutomorphismFunction:

    def __init__(self,function:list,group:CayleyGroup):
        f = []
        for i in function:
            if type(i)==Automorphism:
                f.append(i.phi)
            elif type(i)==tuple:
                if not group.is_automorphism(i):
                    raise ValueError("The transformation is not a valid automorfism for this group.")
                f.append(i)
            else:
                raise TypeError("The argument must be a list of a Automorphisms or tuples.")
            
        self.phi = f
        self.group = group

    def __str__(self):
        return str(self.phi)


#################### HIDDEN COMMANDS ####################

def _obtener_nombre(var_obj):
    for nombre, valor in globals().items():
        if valor is var_obj:
            return nombre
    return None

def _is_closed(tabla):
    n = len(tabla)
    E = _get_elements(tabla)
    if len(E)!=n:
        return False
    for fila in tabla:
        if len(fila) != n:
            return False
    return True

def _check_cayley_group(tabla):
    n = len(tabla)

    # 1. Clausura
    for fila in tabla:
        if len(fila) != n:
            return (False,"ClosingError")
        for x in fila:
            if not (0 <= x < n):
                return (False,"ClosingError")

    # 2. Cada fila y columna debe ser una permutación
    conjunto = set(range(n))

    for fila in tabla:
        if set(fila) != conjunto:
            return (False,f"UniquenessError – row {tabla.index(fila)}")

    for j in range(n):
        columna = {tabla[i][j] for i in range(n)}
        if columna != conjunto:
            return (False,f"UniquenessError – column {j}")

    # 3. Buscar neutro
    neutro = None

    for e in range(n):
        ok = True

        for a in range(n):
            if tabla[e][a] != a:
                ok = False
                break
            if tabla[a][e] != a:
                ok = False
                break

        if ok:
            neutro = e
            break

    if neutro is None:
        return (False,"IndentityError – no identity element found")

    # 4. Inversos
    for a in range(n):
        existe = False

        for b in range(n):
            if tabla[a][b] == neutro and tabla[b][a] == neutro:
                existe = True
                break

        if not existe:
            return (False,f"InverseError – ({a},{b})")

    # 5. Asociatividad
    for a in range(n):
        for b in range(n):
            for c in range(n):

                izquierda = tabla[tabla[a][b]][c]
                derecha = tabla[a][tabla[b][c]]

                if izquierda != derecha:
                    return (False,f"AssociativityError – ({a},{b},{c})")

    return (True,"G is a group")

def _check_explicit_group(elements: list, operation) -> tuple[bool, str]:
    n = len(elements)
    if n == 0:
        return False, "EmptySetError – Group cannot be empty"

    # Mapeo rápido elemento -> índice en O(1) para evitar .index()
    elem_to_idx = {elem: i for i, elem in enumerate(elements)}
    
    if len(elem_to_idx) != n:
        return False, "DuplicateElementError – Elements must be unique"

    # 1. CLAUSURA Y TABULACIÓN RÁPIDA O(n²)
    # Guardamos los resultados como índices para acelerar la comprobación de asociatividad
    cayley_indices = [[0] * n for _ in range(n)]
    
    for i in range(n):
        a = elements[i]
        for j in range(n):
            res = operation(a, elements[j])
            if res not in elem_to_idx:
                return False, f"ClosureError – Result of operation({a}, {elements[j]}) is not in group"
            cayley_indices[i][j] = elem_to_idx[res]

    # 2. ELEMENTO NEUTRO O(n)
    identity_idx = None
    for e_idx in range(n):
        is_identity = True
        for a_idx in range(n):
            # e * a == a y a * e == a
            if cayley_indices[e_idx][a_idx] != a_idx or cayley_indices[a_idx][e_idx] != a_idx:
                is_identity = False
                break
        if is_identity:
            identity_idx = e_idx
            break

    if identity_idx is None:
        return False, "IdentityError – No neutral element found"

    # 3. ELEMENTOS INVERSOS O(n²)
    for a_idx in range(n):
        has_inverse = False
        for b_idx in range(n):
            if cayley_indices[a_idx][b_idx] == identity_idx and cayley_indices[b_idx][a_idx] == identity_idx:
                has_inverse = True
                break
        if not has_inverse:
            return False, f"InverseError – Element '{elements[a_idx]}' has no inverse"

    # 4. ASOCIATIVIDAD O(n³)
    # Se evalúa al final porque es el paso más costoso.
    for a in range(n):
        for b in range(n):
            ab = cayley_indices[a][b]
            for c in range(n):
                # (a * b) * c == a * (b * c)
                if cayley_indices[ab][c] != cayley_indices[a][cayley_indices[b][c]]:
                    return False, f"AssociativityError – Failed for elements ({elements[a]}, {elements[b]}, {elements[c]})"

    return True, "Valid Group"

def _identity(G):
    neutro = None
    for e in range(len(G)):
        ok = True
        for a in range(len(G)):
            if G[e][a] != a:
                ok = False
                break
            if G[a][e] != a:
                ok = False
                break
        if ok:
            neutro = e
            break
    if neutro is None:
        return (False,"IndentityError – no identity element found")
    return neutro

def _subset(G, elements):
    return [[G[r][c] for c in elements] for r in elements]

def _form_subgroup(G,elements):
    # Por propiedades de los subgrupos, solo hace falta comprobar la clausura, el resto se heredan del grupo principal.
    return _is_closed(_reset_renaming(_subset(G,elements)))

def _sign(p):
    """
    Devuelve:
        1  -> permutación par
       -1  -> permutación impar
    """
    inv = 0
    n = len(p)

    for i in range(n):
        for j in range(i+1, n):
            if p[i] > p[j]:
                inv += 1

    return 1 if inv % 2 == 0 else -1

def _renamed_elements(G,L):
    renamed_G = []
    for I in G:
        H = []
        for i in I:
            H.append(L[i])
        renamed_G.append(H)
    return renamed_G

def _get_elements(G):
    """
    Devuelve una lista con todos los elementos/nombres únicos que aparecen
    dentro de la tabla de Cayley G, conservando el orden de primera aparición.
    
    Parámetros:
        G: Lista de listas que representa la tabla de Cayley.
        
    Devuelve:
        Lista con los elementos únicos del grupo.
    """
    # dict.fromkeys() elimina duplicados preservando el orden de aparición
    return list(dict.fromkeys(e for fila in G for e in fila))

def _min_div(n):
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return i
    return n

def _graded_power_set_with_id(G):

    # Principios
    n = len(G)
    m = _min_div(n)
    e = _identity(G)

    def compute_orders(G):

        if not _check_cayley_group(G)[0]:
            raise TypeError

        O = []
        for i in range(len(G)):
            ik = G[i][i]
            o = 1
            while i!=ik:
                ik = G[ik][i]
                o+=1
            O.append(o)

        return O

    orders = compute_orders(G)

    # Posibles órdenes de subgrupos (Lagrange)
    R = []
    for r in range(int(n/m)):
        if n%(r+1)==0:
            R.append(r+1)
    #R.append(n)

    Rm = []
    for r in R[1:]:
        Rm.append(r-1)

    P = []

    # Bucle principal
    for rm in Rm:
        E = []
        for el in range(1,n):
            if (rm+1)%orders[el]==0:
                E.append(el)

        for c in list(combinations(E, rm)):
            P.append(sorted([e]+list(c)))

    return P

def _is_closed_subset(G_cayley, indices):
    """Comprueba clausura en O(|H|^2) usando un conjunto de índices."""
    indices_set = set(indices)
    for i in indices:
        for j in indices:
            if G_cayley[i][j] not in indices_set:
                return False
    return True

def _operate_cosets(G,C,A,B):
    a = A[0]
    b = B[0]
    g = G[a][b]
    for i in C:
        if g in i:
            r = i
            break
    return r

def _order_preserving_permutations(orders):

    classes = {}

    for i, order in enumerate(orders):
        classes.setdefault(order, []).append(i)

    classes = list(classes.values())

    for perms in product(*(permutations(c) for c in classes)):
        result = list(range(len(orders)))

        for domain, image in zip(classes, perms):
            for x, y in zip(domain, image):
                result[x] = y

        yield tuple(result)


#################### PRODUCTS ####################

def direct_product(A:CayleyGroup,B:CayleyGroup):
    G = []
    n, m = A.order(), B.order()

    for a1 in range(n):
        for b1 in range(m):
            g = []
            for a2 in range(n):
                for b2 in range(m):
                    x2 = (a1,b1)
                    g.append( (A.cayley[a1][a2])*m+(B.cayley[b1][b2]) )
            G.append(g)

    names = []
    for a in range(n):
        for b in range(m):
            names.append(str(A.elements[a])+","+str(B.elements[b]))

    return CayleyGroup(G,names,_skip_validation=True)

def direct_power(G:CayleyGroup,n):
    if n<2:
        raise ValueError("n must be at least 2")
    H = G
    for i in range(n-1):
        H = direct_product(H,G)
    return H

def semidirect_product(A: CayleyGroup, B: CayleyGroup, f: AutomorphismFunction):
    G = []
    n, m = A.order(), B.order()

    phi_list = f.phi if hasattr(f, "phi") else f

    for a1 in range(n):
        for b1 in range(m):
            g = []
            phi_b1 = phi_list[b1]
            for a2 in range(n):
                for b2 in range(m):
                    a2_trans = phi_b1[a2]
                    a_prod = A.cayley[a1][a2_trans]
                    b_prod = B.cayley[b1][b2]
                    g.append(a_prod * m + b_prod)
            G.append(g)

    names = []
    for a in range(n):
        for b in range(m):
            names.append(str(A.elements[a]) + "," + str(B.elements[b]))

    return CayleyGroup(G, names, _skip_validation=True)


#################### GENERATORS ####################

def cyclic_group(n):
    elements = range(n)
    G = []
    for i in elements:
        g = []
        for j in elements:
            g.append((i+j)%n)
        G.append(g)
    return CayleyGroup(G,_renaming_C(n),_skip_validation=True)

def symmetric_group(n):

    # Lista de todas las permutaciones
    perms = list(permutations(range(n)))

    # Diccionario permutación -> índice
    index = {p: i for i, p in enumerate(perms)}

    # Tabla
    G = []

    for p in perms:
        fila = []
        for q in perms:
            fila.append(index[tuple(p[i] for i in q)])
        G.append(fila)

    return CayleyGroup(G,_renaming_S(n),_skip_validation=True)

def alternating_group(n):

    perms = [p for p in permutations(range(n)) if _sign(p) == 1]

    index = {p: i for i, p in enumerate(perms)}

    G = []

    for p in perms:
        fila = []
        for q in perms:
            fila.append(index[tuple(p[i] for i in q)])
        G.append(fila)

    return CayleyGroup(G,_renaming_A(n),_skip_validation=True)

def dihedric_group(n):

    elems = [(k,0) for k in range(n)] + [(k,1) for k in range(n)]

    index = {g:i for i,g in enumerate(elems)}

    G = []

    for (k,a) in elems:
        fila = []

        for (l,b) in elems:

            if a == 0:
                m = (k + l) % n
            else:
                m = (k - l) % n

            fila.append(index[(m, a ^ b)])

        G.append(fila)

    return CayleyGroup(G,_renaming_D(n),_skip_validation=True)

def quaternion_group():

    # Multiplicación para la parte positiva
    base = {
        (0,0):( 1,0),
        (0,1):( 1,1),
        (0,2):( 1,2),
        (0,3):( 1,3),

        (1,0):( 1,1),
        (2,0):( 1,2),
        (3,0):( 1,3),

        (1,1):(-1,0),
        (2,2):(-1,0),
        (3,3):(-1,0),

        (1,2):( 1,3),
        (2,3):( 1,1),
        (3,1):( 1,2),

        (2,1):(-1,3),
        (3,2):(-1,1),
        (1,3):(-1,2),
    }

    def multiply(a, b):
        sa, xa = a
        sb, xb = b

        s, x = base[(xa, xb)]

        return (sa * sb * s, x)

    elements = [
        ( 1,0), (-1,0),
        ( 1,1), (-1,1),
        ( 1,2), (-1,2),
        ( 1,3), (-1,3)
    ]

    index = {g: i for i, g in enumerate(elements)}

    G = []

    for a in elements:
        fila = []

        for b in elements:
            fila.append(index[multiply(a, b)])

        G.append(fila)

    return CayleyGroup(G,_renaming_Q8(),_skip_validation=True)

def units_group(n):

    if n <= 1:
        raise ValueError("n debe ser mayor que 1")
    
    elems = [x for x in range(1, n) if gcd(x, n) == 1]
    index = {g: i for i, g in enumerate(elems)}
    
    G = []
    for a in elems:
        fila = []
        for b in elems:
            fila.append(index[(a * b) % n])
        G.append(fila)
        
    return CayleyGroup(G,_renaming_U(n),_skip_validation=True)

def dicyclic_group(n):
    """
    Genera la tabla de Cayley del grupo dicíclico Dic_n (de orden 4n).
    Presentación: <a, x | a^(2n) = 1, x^2 = a^n, x^-1 a x = a^-1>
    """
    if n < 1:
        raise ValueError("n debe ser un entero positivo")
        
    # Elementos representados como pares (k, e) donde 0 <= k < 2n y e en {0, 1}
    # (k, 0) -> a^k
    # (k, 1) -> a^k * x
    elems = [(k, 0) for k in range(2 * n)] + [(k, 1) for k in range(2 * n)]
    index = {g: i for i, g in enumerate(elems)}
    
    G = []
    for (k, a) in elems:
        fila = []
        for (l, b) in elems:
            if a == 0:
                m = (k + l) % (2 * n)
                c = b
            else:
                if b == 0:
                    m = (k - l) % (2 * n)
                    c = 1
                else:
                    m = (k - l + n) % (2 * n)
                    c = 0
            fila.append(index[(m, c)])
        G.append(fila)
        
    return CayleyGroup(G,_renaming_Dic(n),_skip_validation=True)

def tetrahedral_group():
    """
    Full symmetry group of the tetrahedron. Isomorphic to S₄.
    """
    return CayleyGroup(symmetric_group(4).cayley,_skip_validation=True)

def octahedral_group():
    """
    Full symmetry group of the cube/octahedron. Isomorphic to S₄ × C₂.
    """
    return CayleyGroup((symmetric_group(4) * cyclic_group(2)).cayley,_skip_validation=True)

def icosahedral_group():
    """
    Full symmetry group of the icosahedron/dodecahedron. Isomorphic to A₅ × C₂.
    """
    return CayleyGroup((alternating_group(5) * cyclic_group(2)).cayley,_skip_validation=True)


#################### VISUALIZATION AND RENAMING HELPERS ####################

def cayley_table(G, title="", colormap=rainbow, names="", renaming=[]):
    if title=="":
        if _obtener_nombre(G)==None:
            title = f"Cayley table"
        else:
            title = f"Cayley table of {_obtener_nombre(G)}"
    if names=="":
        if len(G)<=20:
            names=True
        else:
            names=False
    if renaming==[]:
        renaming = [rf"${e}$" for e in range(len(G))]

    elements = range(len(G))
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(G, cmap=colormap)
    ax.set_xticks(np.arange(len(elements)))
    ax.set_yticks(np.arange(len(elements)))
    ax.set_xticklabels(renaming)
    ax.set_yticklabels(renaming)

    ax.xaxis.tick_top()
    
    if names:
        for i in range(len(elements)):
            for j in range(len(elements)):
                color_texto = "black"
                ax.text(
                    i,
                    j,
                    f"{renaming[G[j][i]]}",
                    ha="center",
                    va="center",
                    color=color_texto,
                    fontsize=12,
                )
    
    plt.title(title)
    plt.tight_layout()
    plt.show()

def _renaming_C(n):
    r = []
    if n >= 1:
        r.append(r"e")
        if n > 1:
            r.append(r"r")
            if n > 2:
                for i in range(2, n):
                    r.append(rf"r^{{{i}}}")
    return r

def _renaming_S(n):
    """
    Renombrado para S_n en notación de ciclos.
    Coincide con el orden de elementos de generate_group_S(n).
    """
    def tuple_to_cycle(p):
        visited = [False] * len(p)
        cycles = []
        for i in range(len(p)):
            if not visited[i]:
                curr = i
                cycle = []
                while not visited[curr]:
                    visited[curr] = True
                    cycle.append(curr + 1)  # Usamos representación 1-based (1..n)
                    curr = p[curr]
                if len(cycle) > 1:
                    cycles.append("(" + "".join(map(str, cycle)) + ")")
        return "".join(cycles) if cycles else "e"

    perms = list(permutations(range(n)))
    return [tuple_to_cycle(p) for p in perms]

def _renaming_A(n):
    """
    Renombrado para A_n en notación de ciclos.
    Coincide con el orden de elementos de generate_group_A(n).
    """
    def tuple_to_cycle(p):
        visited = [False] * len(p)
        cycles = []
        for i in range(len(p)):
            if not visited[i]:
                curr = i
                cycle = []
                while not visited[curr]:
                    visited[curr] = True
                    cycle.append(curr + 1)
                    curr = p[curr]
                if len(cycle) > 1:
                    cycles.append("(" + "".join(map(str, cycle)) + ")")
        return "".join(cycles) if cycles else "e"

    perms = [p for p in permutations(range(n)) if _sign(p) == 1]
    return [tuple_to_cycle(p) for p in perms]

def _renaming_D(n):
    r = []
    # Rotaciones
    for k in range(n):
        if k == 0:
            r.append(r"e")
        elif k == 1:
            r.append(r"r")
        else:
            r.append(rf"r^{{{k}}}")
            
    # Reflexiones
    for k in range(n):
        if k == 0:
            r.append(r"s")
        elif k == 1:
            r.append(r"rs")
        else:
            r.append(rf"r^{{{k}}}s")
            
    return r

def _renaming_Q8():
    return ["1","-1","i","-i","j","-j","k","-k"]

def _renaming_U(n):
    return [str(x) for x in range(1, n) if gcd(x, n) == 1]

def _renaming_Dic(n):
    names = []
    for k in range(2 * n):
        if k == 0:
            names.append(r"e")
        elif k == 1:
            names.append(r"a")
        else:
            names.append(rf"a^{{{k}}}")
    for k in range(2 * n):
        if k == 0:
            names.append(r"x")
        elif k == 1:
            names.append(r"ax")
        else:
            names.append(rf"a^{{{k}}}x")
    return names

def _reset_renaming(G):
    l = _get_elements(G)
    D = {}
    for i in range(len(l)):
        D[l[i]]=i
    return _renamed_elements(G,D)

