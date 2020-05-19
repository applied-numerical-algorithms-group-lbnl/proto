from st.expr import Index, ConstRef, If, IntLiteral
from st.grid import Grid

# Declare indices
i = Index(0)
j = Index(1)
k = Index(2)

# Declare grid
input = Grid("src", 3)
output = Grid("dest", 3)
coefs = []
for c in range(7):
    coefs.append(ConstRef("coefs[{}]".format(c)))
zero = ConstRef("initToZero")

# Express computation
# output[i, j] is assumed
calc = coefs[0] * input(i, j, k) + \
       coefs[1] * input(i + 1, j, k) + \
       coefs[2] * input(i - 1, j, k) + \
       coefs[3] * input(i, j + 1, k) + \
       coefs[4] * input(i, j - 1, k) + \
       coefs[5] * input(i, j, k + 1) + \
       coefs[6] * input(i, j, k - 1) + \
       output(i, j).assign(If(zero, IntLiteral(0), output(i, j, k)) + calc)

STENCIL = [output]
