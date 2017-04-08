def binom(n,k):
    # homemade binomial coeff computer. Returns float.
    prod = 1.0
    for i in xrange(1,k+1):
            prod *= (n+1-i)
            prod /= i
    return prod

for i in xrange(2,10):print int(binom(2*i,i)),',',
print

from random import randint
from TestAllPolys import GenerateAllMonomials,convertSupportToPoly
def randomPolySupport(DegreeBound, Dimension):
    numberMonoms = int(binom(Dimension + DegreeBound, Dimension))
    powerSetIndex = bin(randint(0,int(2**numberMonoms)))[2:]
    powerSetIndex = powerSetIndex.rjust(numberMonoms,'0')
    monoms = GenerateAllMonomials(DegreeBound, Dimension)
    toReturn = []
    ind = 0
    for monom in monoms:
        if powerSetIndex[ind]=='1':toReturn.append(monom)
        ind+=1
    return toReturn

deg = 3
nvars = 2
print convertSupportToPoly(randomPolySupport(deg, nvars))
