from sympy import symbols

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
