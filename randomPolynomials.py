def binom(n,k):
    # homemade binomial coeff computer. Returns float.
    prod = 1.0
    for i in xrange(1,k+1):
            prod *= (n+1-i)
            prod /= i
    return prod

from random import randint,random
import random
from TestAllPolys import GenerateAllMonomials,convertSupportToPoly
#-------------------------------------------------------------------------------
def randomPolySupports(DegreeBound, Dimension, numPolys=1, weight=0.5):
    """
    Generates random polynomials (support) in Dimension variables and total
    degree <= DegreeBound.
    -- Can get several random polys at a time.
    -- Weight is either the integer number of monomials in each polynomial
       returned, or a float percentage of possible monomials that will occur
       in each polynomial returned.
    """
    numberMonoms = int(binom(Dimension + DegreeBound, Dimension))
    monomsToInclude = 0
    def getPowerSetIndex():
        if type(weight)==float:monomsToInclude = int(round(weight*numberMonoms))
        elif type(weight)==int:monomsToInclude = weight
        else:raise Exception("weight must be of type float or int")
        if monomsToInclude<2:
            #monomsToInclude = 2
            raise Exception("Uh-oh. weight is so small it implies there should \
                  be fewer than 2 monomials.")
        psStr = ['1' for i in xrange(monomsToInclude)]
        psStr += ['0' for i in xrange(numberMonoms-monomsToInclude)]
        random.shuffle(psStr)
        return ''.join(psStr)

    PSindices = [getPowerSetIndex() for i in xrange(numPolys)]
    monoms = GenerateAllMonomials(DegreeBound, Dimension)
    toReturn = [[] for i in xrange(numPolys)]
    index = 0
    for monom in monoms:
        for polyNum in xrange(numPolys):
            if PSindices[polyNum][index]=='1':toReturn[polyNum].append(monom)
        index+=1
    # if we just wanted one, return it,
    # else return all of them   ##nah
    # if len(toReturn)==1:return toReturn[0] 
    return toReturn

#-------------------------------------------------------------------------------
def randomSmallPolys(DegreeBound, Dimension, numMonoms, numPolys=1):
    """
    Generates random polynomials (support) in Dimension variables and total
    degree <= DegreeBound with numMonoms monomials.
    """
    totalNumberMonoms = int(binom(Dimension + DegreeBound, Dimension))
    def randomMonom():
        while True:
            toReturn = tuple([randint(0,DegreeBound+1) for i in xrange(Dimension)])
            if sum(toReturn)<=DegreeBound:return toReturn
    polys = []
    for i in xrange(numPolys):
        thisPoly = {}
        while len(thisPoly)<numMonoms:
            thisPoly[randomMonom()] = 1
        polys.append(thisPoly.keys())
    return polys

if __name__=='__main__':
    deg = 3
    nvars = 2
    print convertSupportToPoly(randomPolySupports(deg, nvars))
    print map(convertSupportToPoly, randomPolySupports(deg, nvars, 3))

