from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout


nautyPath = '/Applications/sage/local/bin/dreadnaut'
with suppress_stdout():
    from phcpy.solver import number_of_symbols, mixed_volume
    from phcpy.polytopes import support
    from phcpy.solver import mixed_volume,linear_product_root_count

#-------------------------------------------------------------------------------
def MakeConstIntoVar(Polys):
    """
    This is a helper function to add the constant term into each monomial list.
    Converting this makes the parsing easier later on.
    """
    for Poly in Polys:
        for Monomial in Poly:
            Monomial.insert(0,0)
            IsConstant = True
            for Term in Monomial:
                if Term != 0:
                    IsConstant = False
                    break
            if IsConstant:
                Monomial[0] = 1
    return Polys

#-------------------------------------------------------------------------------
def CreateCyclicLists(n):
    """
    Create lists of lists of lists in the desired format that represent the
    cyclic n root polynomial system for the given n.
    """
    system = []
    for i in range(n-1):
        equation = []
        for j in range(n):
            mon = [0 for x in range(n)]
            for k in range(i+1):
                mon[(j+k)%n]=1
            equation.append(mon)
        system.append(equation)
    mon1 = [0 for x in range(n)]
    mon2 = [1 for x in range(n)]
    system.append([mon1,mon2])
    return system

#-------------------------------------------------------------------------------
def CreateNautyString(Polys):
    """
    This function takes in the polynomials in the list of lists of lists format
    and converts them into a string that Nauty can read.
    """
    def NautyString(Polys, RootPart, VarPart, MonPart, PolyPart, PowerPart):
    
        def PartString(Partition):
            return str(Partition).replace(',','')[1:-1]
            
        ReturnString = 'c -a -m n=' + str(len(Polys)) + ' g '
        for Poly in Polys:
            for Term in Poly:
                ReturnString += str(Term) + ' '
        
            ReturnString += ';'
        ReturnString = ReturnString[:-1] + '. f = [ 0 | ' + str(RootPart) + ' | ' 
        ReturnString += PartString(VarPart) + ' | ' + PartString(MonPart) 
        ReturnString += ' | ' + PartString(PolyPart) + ' | '
        for Part in PowerPart:
            ReturnString += PartString(Part) + ' | '
        ReturnString +=  '] x @ b'
        return ReturnString
    
    SystemAsLists = []
    n = len(Polys[0][0])

    for i in range(n):
        SystemAsLists.append([])
        
    Monomials = []
    Polynomials = []
    Variables = list(range(n))
    Variables.remove(0)
    NewNodeRef = n
    TermToNode = {}
    SystemNode = NewNodeRef
    SystemAsLists.append([])
    NewNodeRef += 1
    for i in range(len(Polys)):
        PolyNode = NewNodeRef
        Polynomials.append(NewNodeRef)
        SystemAsLists.append([SystemNode])
        NewNodeRef += 1
        for j in range(len(Polys[i])):
            MonomialNode = NewNodeRef
            Monomials.append(NewNodeRef)
            SystemAsLists[PolyNode].append(MonomialNode)
            SystemAsLists.append([])
            NewNodeRef += 1
            for k in range(len(Polys[i][j])):
                if Polys[i][j][k] != 0:
                    if (k, Polys[i][j][k]) not in TermToNode:
                        TermToNode[(k, Polys[i][j][k])] = NewNodeRef
                        SystemAsLists.append([])
                        NewNodeRef += 1
                        SystemAsLists[k].append(TermToNode[(k, Polys[i][j][k])])
                    SystemAsLists[MonomialNode].append(TermToNode[(k, Polys[i][j][k])])
    #print "This is the system converted into our pre-graph list format:"
    #print SystemAsLists
    
    PartList = [[]]
    ExponentsToPartition = list(TermToNode.keys())
    
    ExponentsToPartition.sort(key=lambda x: x[1])
    CurrentPower = ExponentsToPartition[0][1]
    global PowerString
    PowerString = str(CurrentPower) + ' '
    for i in range(len(ExponentsToPartition)):
        PossiblyNewPower = ExponentsToPartition[i][1]
        if CurrentPower != PossiblyNewPower:
            CurrentPower = PossiblyNewPower
            PowerString += str(CurrentPower) + ' '
    if len(ExponentsToPartition) > 0:
        MinValue = ExponentsToPartition[0][1]
    for Key in ExponentsToPartition:
        if Key[1] != MinValue:
            PartList.append([])
            MinValue = Key[1]
        PartList[-1].append(TermToNode[Key])
    return NautyString(SystemAsLists, SystemNode, Variables, Monomials, Polynomials, PartList)

#-------------------------------------------------------------------------------
def GetNautyOutput(CanonizedLists):
    """
    This function takes a polynomial in the desired format, and then it performs
    the nauty call and parses the output to create the unique string.
    """
    from subprocess import Popen, PIPE, STDOUT
    p = Popen([nautyPath], stdout = PIPE, stdin = PIPE, stderr = STDOUT, universal_newlines = True)
    Output = p.communicate(input = CanonizedLists)[0]
    return Output

#-------------------------------------------------------------------------------
def GetUniqueString(Output):
    return Output[Output.find(':') - 4:] + PowerString

#-------------------------------------------------------------------------------
def canonize(pols):
    dim = number_of_symbols(pols)
    
    supp = [[list(mon) for mon in support(dim,pol)] for pol in pols]
    nautyString = GetNautyOutput(CreateNautyString(MakeConstIntoVar(supp)))
    groupSize = GetGroupSize(nautyString,supp)
    uniqueString = GetUniqueString(nautyString)
    return uniqueString,groupSize

#-------------------------------------------------------------------------------
def GetGroupSize(nautyString,supp):
    # finds the system's symmetry group size (after dividing out by the 
    # "same polynomial" symmetry)
    from math import factorial
    grpStart = nautyString.find('grpsize=') + 8 
    grpEnd = nautyString.find(';',grpStart) 
    groupSize = int(nautyString[grpStart:grpEnd])
    polySymmetry = reduce(lambda a,b:a*b,map(factorial,SystemSupportPartition(supp)))
    return groupSize/polySymmetry

#-------------------------------------------------------------------------------
def SystemSupportPartition(supp):
    count = {}
    polys = [tuple([tuple(m) for m in p]) for p in supp]
    for poly in polys:
        if poly in count:count[poly] += 1
        else:count[poly]=1
    return count.values()

#-------------------------------------------------------------------------------
def GenerateAllMonomials(DegreeBound, Dimension):
    from itertools import ifilter,product
    degBoundFn = lambda a:sum(a)<=DegreeBound
    return ifilter(degBoundFn,product(range(DegreeBound+1),repeat=Dimension))

#-------------------------------------------------------------------------------
def GenerateAllSystems(DegreeBound, Dimension, includeConst=False, includeMonoms=False):
    """
    Generate all systems in Dimension variables with degree <= DegreeBound.
    includeConst determines whether we include the constant polynomial.
    includeMonoms determines whether we include single-monomial polynomials.
    includeMonoms being false overrides includeConst, as constants are special
    cases of monomials.
    Returns an iterator, rather than taking up tons of memory (I hope).
    """
    from itertools import chain, combinations,combinations_with_replacement
    def binom(n,k):
        # homemade binomial coeff computer
        prod = 1.0
        for i in xrange(1,k+1):
                prod *= (n+1-i)
                prod /= i
        return prod
    numberOfMonoms = int(binom(DegreeBound+Dimension,Dimension))
    monoms = GenerateAllMonomials(DegreeBound, Dimension)
    if includeMonoms:includer = 1
    elif includeConst:includer = 0
    else:includer = 2
    polys = chain.from_iterable(combinations(monoms,n) for n in xrange(includer,numberOfMonoms))
    systs = combinations_with_replacement(polys,Dimension)
    return systs

#-------------------------------------------------------------------------------
def convertSupportToPoly(supp):
    """
    Takes in a polynomial,
    returns it as a polynomial string for use by phcpy
    """
    nvars = len(supp[0])
    variables = ['x'+str(i) for i in xrange(1,nvars+1)]
    thisPoly = []
    for monomial in supp:
        thisPoly.append('*'.join([variables[i]+'^'+str(monomial[i]) for i in xrange(nvars)]))
    return ' + '.join(thisPoly)+ ';'

if __name__=='__main__':
    mixedVolumes = []
    linearProdCounts = []
    grpSizes = []
    numVars = 2
    degreeBound = 2
    Supports = GenerateAllSystems(degreeBound, numVars)
    for syst in Supports:
        pols = map(convertSupportToPoly, syst)
        mixedVolumes.append(mixed_volume(pols))
        linearProdCounts.append(linear_product_root_count(pols, silent=True))
        grpSizes.append(canonize(pols)[1])

    fileEnding = 'dim='+str(numVars)+'_deg='+str(degreeBound)
    with open('MV_'+fileEnding,'w') as f:
        f.write(str(mixedVolumes))
    with open('linProds_'+fileEnding,'w') as f:
        f.write(str(linearProdCounts))
    with open('groupSizes_'+fileEnding,'w') as f:
        f.write(str(grpSizes))
