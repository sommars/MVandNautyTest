nautyPath = '/Applications/sage/local/bin/dreadnaut'
from phcpy.solver import number_of_symbols, mixed_volume
from phcpy.polytopes import support

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
    print groupSize,polySymmetry,'<---'
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
def MixedVolume(pols):
    return mixed_volume(pols)

#-------------------------------------------------------------------------------
def MonsToPols(Mons,includeConst=False,includeMonoms=False):
    # return an iterable of polynomials
    from itertools import chain, combinations
    if includeConst:includer = 0
    elif includeMonoms:includer = 1
    else:includer = 2
    return chain.from_iterable(combinations(Mons,n) for n in range(includer,len(Mons)+1))
    return list(toReturn)

#-------------------------------------------------------------------------------
def GenerateTwoDimMonsOfDegree(d):
    Pts = []
    for i in xrange(0,d+1):
        for j in xrange(0,(d+1)):
           if i + j <= d:
              Pts.append([i,j])
    Pols = MonsToPols(Pts)
    Result = []
    for Pol in Pols:
        for Pol2 in Pols:
            Result.append([Pol,Pol2])
    return Result
    
#-------------------------------------------------------------------------------
def GenerateThreeDimMonsOfDegree(d):
    Pts = []
    for i in xrange(0,d+1):
        for j in xrange(0,(d+1)):
            for k in xrange(0,(d+1)):
                if i + j + k <= d:
                    Pts.append([i,j,k])
    Pts.sort()
    Pts = Pts[1:]
    Result = []
    for Pt in Pts:
        for Pt2 in Pts:
            for Pt3 in Pt2:
                Result.append([Pt,Pt2, Pt3])
    return Result

#-------------------------------------------------------------------------------
def GenerateFourDimMonsOfDegree(d):
    Pts = []
    for i in xrange(0,d+1):
        for j in xrange(0,(d+1)):
            for k in xrange(0,(d+1)):
                for l in xrange(0,(d+1)):
                    if i + j + k + l<= d:
                        Pts.append([i,j,k,l])
    Pts.sort()
    return Pts

#-------------------------------------------------------------------------------
def GenerateAllSystems(DegreeBound, Dimension, includeConst=False, includeMonoms=False):
    from itertools import product, ifilter, chain, combinations
    from scipy.special import binom
    degBoundFn = lambda a:sum(a)<=DegreeBound
    monoms = list(ifilter(degBoundFn,product(range(DegreeBound+1),repeat=Dimension)))
    if includeConst:includer = 0
    elif includeMonoms:includer = 1
    else:includer = 2
    print int(binom(DegreeBound,Dimension))+2
    polys = chain.from_iterable(combinations(monoms,n) for n in xrange(includer,int(binom(DegreeBound,Dimension))+2))
    return list(polys)

#-------------------------------------------------------------------------------
def convertSupportToPolys(supp):
    polys = []
    nvars = len(supp[0][0])
    variables = ['x'+str(i) for i in xrange(1,nvars+1)]
    for poly in supp:
        thisPoly = []
        for monomial in poly:
            thisPoly.append('*'.join([variables[i]+'^'+str(monomial[i]) for i in xrange(nvars)]))
        polys.append(' + '.join(thisPoly)+ ';')
    return polys

"""
for i in GenerateAllSystems(2,2):print i
#for i in convertSupportToPolys(GenerateAllSystems(2,2)):print i

Supports = GenerateAllSystems(2,2)
print(len(Supports))
#for i in convertSupportToPolys(Supports):print i
print Supports[0]
#for syst in Supports:print convertSupportToPolys(syst)
systs = map(convertSupportToPolys,Supports)
print systs[10]
#print canonize(systs[100])
"""
