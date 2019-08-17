from sympy import symbols

a = "12;1"
print(a.split(';'))

x = symbols("x")
y = symbols("y")


test = (x, y)

print(all([x in test for x in test]))
