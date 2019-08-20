from sympy import symbols
import sympy
a = "12;1"

x = symbols("x")
y = symbols("y")


test = (x, y)



x = False
y = False
z = False

def myt():
    return not(x or y or z)

print(myt())

x= True

print(myt())

eq = sympy.sympify("A**2*k1*k4*k7*k_3 + A*B*k1*k3*k6*k7 + A*B*k1*k3*k6*k_4")
print(sympy.factor(eq, deep=True))