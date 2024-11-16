#########################################################################################################
#                                                                                                       #
#                    TALLER APLICACIONES ÁLGEBRA - JUAN DAVID SIERRA OROZCO                             #
#                                                                                                       #
#                                                                                                       #
#########################################################################################################    

#PUNTO 1

import numpy as np
from fractions import Fraction
from scipy.linalg import null_space

def parse_compound(compound):
    """Analiza un compuesto químico en un diccionario de elementos y sus cantidades."""
    elementos = {}
    i = 0
    while i < len(compound):
        # Manejar paréntesis
        if compound[i] == '(':
            close_paren = compound.find(')', i)
            subcompound = compound[i+1:close_paren]
            i = close_paren + 1
            # Verificar si hay un número después del paréntesis
            num = ''
            while i < len(compound) and compound[i].isdigit():
                num += compound[i]
                i += 1
            multiplier = int(num) if num else 1
            sub_elements = parse_compound(subcompound)
            for element, count in sub_elements.items():
                elementos[element] = elementos.get(element, 0) + count * multiplier
            continue
            
        # Manejar elementos regulares
        if compound[i].isupper():
            elemento = compound[i]
            i += 1
            while i < len(compound) and compound[i].islower():
                elemento += compound[i]
                i += 1
            
            # Obtener número después del elemento
            num = ''
            while i < len(compound) and compound[i].isdigit():
                num += compound[i]
                i += 1
            elementos[elemento] = elementos.get(elemento, 0) + (int(num) if num else 1)
        else:
            i += 1
            
    return elementos

def balance_equation(reactivos, productos):
    """
    Balancea una ecuación química usando álgebra lineal.
    
    Args:
        reactivos: Lista de compuestos reactivos
        productos: Lista de compuestos productos
    
    Returns:
        Lista de coeficientes para la ecuación balanceada
    """
    # Obtener todos los elementos únicos
    elementos = set()
    for compound in reactivos + productos:
        elementos.update(parse_compound(compound).keys())
    elementos = sorted(list(elementos))
    
    # Crear la matriz de coeficientes
    num_compounds = len(reactivos) + len(productos)
    matriz = np.zeros((len(elementos), num_compounds))
    
    # Llenar la matriz con los reactivos (coeficientes positivos)
    for i, compound in enumerate(reactivos):
        parsed = parse_compound(compound)
        for elemento, count in parsed.items():
            row = elementos.index(elemento)
            matriz[row, i] = count
            
    # Llenar la matriz con los productos (coeficientes negativos)
    for i, compound in enumerate(productos):
        parsed = parse_compound(compound)
        for elemento, count in parsed.items():
            row = elementos.index(elemento)
            matriz[row, i + len(reactivos)] = -count
    
    # Imprimir la matriz de coeficientes
    print("Matriz de Coeficientes:")
    print(matriz)
    
    # Encontrar el espacio nulo usando scipy
    nullspace = null_space(matriz)
    
    if nullspace.size == 0:
        return None
    
    # Obtener la primera solución del espacio nulo
    solucion = nullspace[:, 0]
    
    # Escalar la solución a los valores enteros más pequeños
    min_positive = np.min(np.abs(solucion[np.nonzero(solucion)]))
    scaled_solution = solucion / min_positive
    coeficientes = [int(round(x)) for x in scaled_solution]
    
    # Hacer todos los coeficientes positivos multiplicando por -1 si es necesario
    if coeficientes[0] < 0:
        coeficientes = [-x for x in coeficientes]
    
    # Reducir los coeficientes por su MCD
    mcd = np.gcd.reduce(np.abs(coeficientes))
    if mcd > 1:
        coeficientes = [x // mcd for x in coeficientes]
    
    return coeficientes

def format_balanced_equation(reactivos, productos, coeficientes):
    """Formatea la ecuación balanceada como una cadena de texto."""
    len_reactivos = len(reactivos)
    
    # Formatear reactivos
    terminos_reactivos = []
    for i, (coef, compound) in enumerate(zip(coeficientes[:len_reactivos], reactivos)):
        term = f"{coef if coef > 1 else ''}{compound}"
        terminos_reactivos.append(term)
    
    # Formatear productos
    terminos_productos = []
    for i, (coef, compound) in enumerate(zip(coeficientes[len_reactivos:], productos)):
        term = f"{coef if coef > 1 else ''}{compound}"
        terminos_productos.append(term)
    
    return f"{' + '.join(terminos_reactivos)} → {' + '.join(terminos_productos)}"

def solve_equation(equation_str):
    """Resuelve una ecuación química dada como una cadena de texto."""
    # Separar en reactivos y productos
    reactivos_str, productos_str = equation_str.split('→')
    
    # Analizar reactivos y productos
    reactivos = [x.strip() for x in reactivos_str.split('+')]
    productos = [x.strip() for x in productos_str.split('+')]
    
    # Balancear la ecuación
    coeficientes = balance_equation(reactivos, productos)
    
    if coeficientes is None:
        return "No se encontró solución"
    
    # Formatear la ecuación balanceada
    return format_balanced_equation(reactivos, productos, coeficientes)

# Lista de ecuaciones a resolver
ecuaciones = [
    "FeS2 + O2 → Fe2O3 + SO2",
    "CO2 + H2O → C6H12O6 + O2",
    "C4H10 + O2 → CO2 + H2O",
    "C5H11OH + O2 → H2O + CO2",
    "HClO4 + P4O10 → H3PO4 + Cl2O7",
    "Na2CO3 + C + N2 → NaCN + CO",
    "C2H2Cl4 + Ca(OH)2 → C2HCl3 + CaCl2 + H2O"
]

# Resolver cada ecuación
print("Ecuaciones Químicas Balanceadas:")
print("-" * 50)
for i, equation in enumerate(ecuaciones, 1):
    print(f"{i}. {equation}")
    balanced = solve_equation(equation)
    print(f"Balanceada: {balanced}")
    print()
