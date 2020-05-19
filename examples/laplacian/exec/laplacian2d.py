from st.expr import Index, ConstRef, If, IntLiteral
from st.grid import Grid

# Declare indices
i = Index(0)
j = Index(1)

# Declare grid
input = Grid("src", 2)
output = Grid("dest", 2)
coefs = []
for c in range(5):
    coefs.append(ConstRef("coefs[{}]".format(c)))
zero = ConstRef("initToZero")

# Express computation
# output[i, j] is assumed
calc = coefs[0] * input(i, j) + \
       coefs[1] * input(i + 1, j) + \
       coefs[2] * input(i - 1, j) + \
       coefs[3] * input(i, j + 1) + \
       coefs[4] * input(i, j - 1)
output(i, j).assign(If(zero, IntLiteral(0), output(i, j)) + calc)

STENCIL = [output]
