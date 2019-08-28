from sympy import symbols
import sympy




# k_ib = symbols("K_iB")
# k_ma = symbols("K_mA")
# k_ma2 = symbols("K_mA2")
# k_mb = symbols("K_mB")
# k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8 =symbols(["k_"+str((x+1)*-1) for x in range(8)])
# k1, k2, k3, k4, k5, k6, k7, k8 = symbols("k1:9")
# #k1, (k_1 + k2)/k_ma)

# eq = sympy.Eq((k_1 + k2)/k_ma, k1)
# #sympy.pprint(eq)
# #sympy.pprint(sympy.solve(eq, k_ma))

from pylatex.base_classes import Environment
from pylatex.package import Package
from pylatex import Document, Section
from pylatex.utils import NoEscape
# ns = dict()

# ns["k_-1"] = symbols("k_-1")
# sympy.pprint(ns["k_-1"])
# sympy.pprint(sympy.sympify("(k_-1 + k2)/k_ma", locals=ns))
# sympy.pprint(sympy.sympify("(k_-1 + k2)/k_ma"))

txt = f'{sympy.latex(self.from_state)} & \\xrightarrow{{{sympy.latex(self.rate)}}} & {sympy.latex(self.from_state)}'
print(txt)