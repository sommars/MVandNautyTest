from randomPolynomials import randomPolySupports,binom
from TestAllPolys import convertSupportToPoly,canonize
from phcpy.solver import mixed_volume as MV
from phcpy.solver import linear_product_root_count as LPRC
def bezout(syst):
    # for a poly syst as list of supports, return total deg
    degrees = map(lambda p:max(map(sum,p)),syst)
    return reduce(lambda a,b:a*b,degrees)

if __name__=='__main__':
#def asdf():
    numVars = 4
    degreeBound = 16
    #weightVal = 0.151
    weightVal = 21
    numberTrials = 400
    print "deg:",degreeBound," #vars:",numVars
    print "poly density:",weightVal,"#monomials:",weightVal#int(round(weightVal*binom(degreeBound+numVars,numVars)))
    print "#trials:",numberTrials
    data = []
    for i in xrange(numberTrials):
        syst = randomPolySupports(degreeBound, numVars, numPolys=numVars, weight=weightVal)
        polys = map(convertSupportToPoly,syst)
        #newDatum = str((bezout(syst), LPRC(polys, silent=True), MV(polys), canonize(polys)[1]))
        newDatum = MV(polys)*1.0/bezout(syst)
        data.append(newDatum)
    #print ', '.join(data)
    print "Average ratio of MV:Bezout:",sum(data)/len(data)
    fileEnding = 'dim='+str(numVars)+'_deg='+str(degreeBound)+'_random'
    with open('data_'+fileEnding,'w') as f:
        f.write('bezout,LPRC,MV,gpSize\n')
        f.write('\n'.join(data))

#import cProfile
#cProfile.run('asdf()')
