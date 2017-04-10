def binom(n,k):
    # homemade binomial coeff computer. Returns float.
    prod = 1.0
    for i in xrange(1,k+1):
            prod *= (n+1-i)
            prod /= i
    return prod

from random import randint,random
from TestAllPolys import GenerateAllMonomials,convertSupportToPoly
'''
#-------------------------------------------------------------------------------
def randomPolySupports(DegreeBound, Dimension, numPolys=1):
    """
    Generates random polynomials in Dimension variables and total
    degree <= DegreeBound. Equivalent to finding a random element of
    the power set of all monomials, so if there are M of them, we pick
    a random int 0<=x<=M and use its binary representation to choose
    from the list of monomials. Can get several random polys at a time.
    Weight is the probability that a monomial gets chosen.
    """
    numberMonoms = int(binom(Dimension + DegreeBound, Dimension))
    def getPowerSetIndex():
        powerSetIndex = bin(randint(0,int(2**numberMonoms)))[2:]
        return powerSetIndex.rjust(numberMonoms,'0')
    PSindices = [getPowerSetIndex() for i in xrange(numPolys)]
    monoms = GenerateAllMonomials(DegreeBound, Dimension)
    toReturn = [[] for i in xrange(numPolys)]
    index = 0
    for monom in monoms:
        for polyNum in xrange(numPolys):
            if PSindices[polyNum][index]=='1':toReturn[polyNum].append(monom)
        index+=1
    # if we just wanted one, return it,
    # else return all of them
    if len(toReturn)==1:return toReturn[0] 
    return toReturn
'''

#-------------------------------------------------------------------------------
def randomPolySupports(DegreeBound, Dimension, numPolys=1, weight=0.5):
    """
    Generates random polynomials in Dimension variables and total
    degree <= DegreeBound. Equivalent to finding a random element of
    the power set of all monomials, so if there are M of them, we pick
    a random int 0<=x<=M and use its binary representation to choose
    from the list of monomials. Can get several random polys at a time.
    Weight is the probability that a monomial gets chosen.
    """
    numberMonoms = int(binom(Dimension + DegreeBound, Dimension))
    monoms = GenerateAllMonomials(DegreeBound, Dimension)
    toReturn = [[] for i in xrange(numPolys)]
    for monom in monoms:
        for polyNum in xrange(numPolys):
            if random()<weight:toReturn[polyNum].append(monom)
    # if we just wanted one, return it,
    # else return all of them
    if len(toReturn)==1:return toReturn[0] 
    return toReturn


if __name__=='__main__':
    deg = 3
    nvars = 2
    print convertSupportToPoly(randomPolySupports(deg, nvars))
    print map(convertSupportToPoly, randomPolySupports(deg, nvars, 3))

