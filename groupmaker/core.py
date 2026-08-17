#————————————————————————————————————————————————————————————#
#                                                            #
#                        GroupMaker                          #
#          A Python library for finite group theory          #
#                                                            #
#                  by: Mario Sultan Romero                   #
#                                                            #
#————————————————————————————————————————————————————————————#


#################### IMPORTS ####################

import  numpy               as      np
import  matplotlib.pyplot   as      plt
import  matplotlib.colors   as      mcolors
from    itertools           import  permutations, chain, combinations, product
from    collections         import  Counter
from    math                import  gcd


#################### DEFINITIONS ####################

tgl_color = "#2A646E"
white = mcolors.LinearSegmentedColormap.from_list("white", ["white", "white"])
tgl = mcolors.LinearSegmentedColormap.from_list("tgl", ["white", tgl_color])
rainbow = mcolors.LinearSegmentedColormap.from_list("rainbow", ["#FF0000","#FF9100","#F2DA00","#2CDB00","#00DAE9","#1869FF"])#,"#7648FF","#D12BFF"])


#################### GROUP CLASS ####################

class Group:

    def __init__(self, cayley, names=None):

        """
        if type(cayley)!=list[list[int]]:
            return ValueError("The group must be a list of lists of integers")
        """

        valid, message = is_group(cayley)

        if not valid:
            raise ValueError(message)

        if names is None:
            names = [i for i in range(len(cayley))]

        if len(names) != len(cayley):
            raise ValueError("The number of names must match the group order")

        self.cayley = cayley
        self.names = names

    def __str__(self):
        G = []
        for i in range(self.order()):
            g = []
            for j in range(self.order()):
                g.append(self.names[self.cayley[i][j]])
            G.append(g)
        s = str(G)
        return f"{G}"

    def _print_group(self):
        print(f"Group({self.cayley},{self.names})")

    def order(self):
        return len(self.cayley)

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

    def proper_subgroups(self):
        G = self.cayley
        E = self.names
        L = []
        n = self.order()
        for i in _graded_power_set_with_id(G):
            ii = []
            for j in i:
                ii.append(E[j])
            if form_subgroup(G,i):
                L.append(Group(_reset_renaming(subset(G,i)),ii))         
        return L

    def subgroups(self):
        G = self.cayley
        S = self.proper_subgroups()[:]
        S.insert(0,Group([[0]]))
        S.append(self)
        return S

    def cayley_table(self, title="", colormap=rainbow, names=None):

        if names==None:
            if len(self.cayley)<=20:
                names=True
            else:
                names=False


        elements = range(len(self.cayley))
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(self.cayley, cmap=colormap)
        ax.set_xticks(np.arange(len(elements)))
        ax.set_yticks(np.arange(len(elements)))
        ax.set_xticklabels(self.names)
        ax.set_yticklabels(self.names)
        ax.xaxis.tick_top()
        
        if names:
            for i in range(len(elements)):
                for j in range(len(elements)):
                    color_texto = "black"
                    ax.text(
                        i,
                        j,
                        f"{self.names[self.cayley[j][i]]}",
                        ha="center",
                        va="center",
                        color=color_texto,
                        fontsize=12,
                    )
        
        plt.title(title)
        plt.tight_layout()
        plt.show()

    def delete_names(self):
        self.names = [i for i in range(len(self.cayley))]




#class Element:

#class Automorphism:
    

#################### HIDDEN COMMANDS ####################

def _obtener_nombre(var_obj):
    for nombre, valor in globals().items():
        if valor is var_obj:
            return nombre
    return None

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

def _renaming_C(n):
    r = []
    if n >= 1:
        r.append(r"$e$")
        if n > 1:
            r.append(r"$r$")
            if n > 2:
                for i in range(2, n):
                    r.append(rf"$r^{{{i}}}$")
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
            r.append(r"$e$")
        elif k == 1:
            r.append(r"$r$")
        else:
            r.append(rf"$r^{{{k}}}$")
            
    # Reflexiones
    for k in range(n):
        if k == 0:
            r.append(r"$s$")
        elif k == 1:
            r.append(r"$rs$")
        else:
            r.append(rf"$r^{{{k}}}s$")
            
    return r

def _renaming_Q8():
    return ["1","-1","i","-i","j","-j","k","-k"]

def _renaming_U(n):
    return [str(x) for x in range(1, n) if gcd(x, n) == 1]

def _renaming_Dic(n):
    names = []
    for k in range(2 * n):
        if k == 0:
            names.append(r"$e$")
        elif k == 1:
            names.append(r"$a$")
        else:
            names.append(rf"$a^{{{k}}}$")
    for k in range(2 * n):
        if k == 0:
            names.append(r"$x$")
        elif k == 1:
            names.append(r"$ax$")
        else:
            names.append(rf"$a^{{{k}}}x$")
    return names

def _reset_renaming(G):
    l = get_elements(G)
    D = {}
    for i in range(len(l)):
        D[l[i]]=i
    return renamed_elements(G,D)

def _min_div(n):
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return i
    return n

def _power_set_with_id(G):
    n = order(G)
    m = _min_div(n)
    lista = get_elements(G)
    e = identity(G)
    lista.remove(e)
    P = list(list(c) for c in chain.from_iterable(combinations(lista, r) for r in range(int(n/m))))
    for i in range(len(P)):
        P[i].insert(0,e)
        #P[i] = sorted(P[i])
    return P

def _graded_power_set_with_id(G):

    # Principios
    n = order(G)
    m = _min_div(n)
    e = identity(G)
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

def _operate_automorphisms(a,b):
    aob = []
    for x in b:
        aob.append(a[x])
    return tuple(aob)

def _enumeration_dict(L):
    D = {}
    for i in range(len(L)):
        D[tuple(L[i])]=i
    return D

def _cosets(G,H):
    C = []
    for g in get_elements(G):
        gH = coset(G,H,g)
        if gH not in C:
            C.append(gH)
    return C

def _operate_cosets(G,C,A,B):
    a = A[0]
    b = B[0]
    g = G[a][b]
    for i in C:
        if g in i:
            r = i
            break
    return r


#################### GENERAL COMMANDS ####################

def renamed_elements(G,L):
    renamed_G = []
    for I in G:
        H = []
        for i in I:
            H.append(L[i])
        renamed_G.append(H)
    return renamed_G

def get_elements(G):
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


#################### BASIC COMMANDS ####################

def subset(G, elements):
    return [[G[r][c] for c in elements] for r in elements]

def identity(G):
    neutro = None
    for e in range(order(G)):
        ok = True
        for a in range(order(G)):
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

def is_closed(tabla):
    n = len(tabla)
    E = get_elements(tabla)
    if len(E)!=n:
        return False
    for fila in tabla:
        if len(fila) != n:
            return False
    return True

def is_group(tabla):
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

def order(G):
    if is_group(G)[0]:
        return len(G)
    else:
        raise ValueError("El argumento introducido no es un grupo")

def compute_orders(G):

    if not is_group(G)[0]:
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

def count_orders(G):
    l = sorted(compute_orders(G))
    return dict(Counter(l))

def is_cyclic(G):
    if is_group(G)[0]:
        return order(G) in compute_orders(G)
    else:
        raise ValueError("El argumento introducido no es un grupo")

def center(G):
    if not is_group(G)[0]:
        raise ValueError("El argumento introducido no es un grupo")
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

def is_abelian(G):
    return order(G)==len(center(G))

def form_subgroup(G,elements):
    # Por propiedades de los subgrupos, solo hace falta comprobar la clausura, el resto se heredan del grupo principal.
    return is_closed(_reset_renaming(subset(G,elements)))

def proper_subgroups(G):
    L = []
    n = order(G)
    for i in _graded_power_set_with_id(G):
        if form_subgroup(G,i):
            L.append(i)   
    #L.append(list(range(n)))       
    return L

def subgroups(G):
    S = proper_subgroups(G)[:]
    S.insert(0,[0])
    S.append(get_elements(G))
    return S

def coset(G,H,element,side="left"):
    L = []
    for h in H:
        if side=="left":
            ah = G[element][h]
        elif side=="right":
            ah = G[h][element]
        if ah not in L:
            L.append(ah)
    return sorted(L)

def is_normal(G,H):
    normal = True
    for g in get_elements(G):
        if coset(G,H,g,side="left")!=coset(G,H,g,side="right"):
            normal=False
            break
    return normal

def normal_subgroups(G):
    NS = []
    for H in subgroups(G):
        if is_normal(G,H):
            NS.append(H)
    return NS

def is_simple(G):
    return normal_subgroups(G)==[[0],get_elements(G)]

def quotient_group(G,H):
    C = _cosets(G,H)
    D = _enumeration_dict(C)
    Q = []
    for i in range(len(C)):
        A = []
        for j in range(len(C)):
            A.append(D[tuple(_operate_cosets(G,C,C[i],C[j]))])
        Q.append(A)
    return (Q,D)

def is_automorphism(G,φ):
    # a es una n-tupla
    n = len(φ)
    if not is_group(G):
        raise TypeError ("The first argument is not a group")
    if type(φ)!=tuple:
        raise TypeError ("The second argument is not a tuple")
    
    # Biyectividad
    for i in range(n):
        if i not in φ:
            return False

    #Elemento neutro
    if φ[0]!=0:
        return False

    #Preservación de la LCI
    for i in range(len(G)):
        for j in range(len(G)):
            if G[φ[i]][φ[j]]!=φ[G[i][j]]:
                return False
    
    return True

def automorphisms(G):
    g = tuple(range(len(G)))
    Auts = []
    for φ in permutations(g):
        if is_automorphism(G,φ):
            Auts.append(φ)
    return Auts

def automorphism_group(G):
    auts = automorphisms(G)
    n = len(auts)
    D = _enumeration_dict(auts)
    Aut = []
    for i in range(n):
        a = []
        for j in range(n):
            a.append(D[_operate_automorphisms(auts[i],auts[j])])
        Aut.append(a)
    return (Aut,D)


#################### GENERATORS ####################

def direct_product(A,B):

    if (not is_group(A)[0]) or (not is_group(B)[0]):
        raise TypeError

    G = []

    n, m = len(A), len(B)

    for a1 in range(n):
        for b1 in range(m):
            g = []
            for a2 in range(n):
                for b2 in range(m):
                    x2 = (a1,b1)
                    g.append( (A[a1][a2])*len(B)+(B[b1][b2]) )
            G.append(g)

    return G

def direct_power(G,n):
    H = [[0]]
    for i in range(n):
        H = direct_product(H,G)
    return H

def cyclic_group(n):
    elements = range(n)
    G = []
    for i in elements:
        g = []
        for j in elements:
            g.append((i+j)%n)
        G.append(g)

    return Group(G,_renaming_C(n))

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

    return Group(G,_renaming_S(n))

def alternating_group(n):

    perms = [p for p in permutations(range(n)) if _sign(p) == 1]

    index = {p: i for i, p in enumerate(perms)}

    G = []

    for p in perms:
        fila = []
        for q in perms:
            fila.append(index[tuple(p[i] for i in q)])
        G.append(fila)

    return Group(G,_renaming_A(n))

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

    return Group(G,_renaming_D(n))

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

    return Group(G,_renaming_Q8())

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
        
    return Group(G,_renaming_U(n))

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
        
    return Group(G,_renaming_Dic(n))


#################### VISUALIZATION AND RENAMING HELPERS ####################

def print_group(G):
    for i in G:
        for j in i:
            print(j,end="\t")
        print()

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

