from randomPolynomials import randomPolySupports
from TestAllPolys import convertSupportToPoly,canonize
from phcpy.solver import mixed_volume as MV
from phcpy.solver import linear_product_root_count as LPRC
def bezout(syst):
    # for a poly syst as list of supports, return total deg
    degrees = map(lambda p:max(map(sum,p)),syst)
    return reduce(lambda a,b:a*b,degrees)

if __name__=='__main__':
#def asdf():
    numVars = 5
    degreeBound = 5
    probOfMonom = 0.3
    data = []
    for i in xrange(24):
        syst = randomPolySupports(degreeBound, numVars, numPolys=numVars, weight=probOfMonom)
        polys = map(convertSupportToPoly,syst)
        newDatum = str((bezout(syst), LPRC(polys, silent=True), MV(polys), canonize(polys)[1]))
        data.append(newDatum)
    fileEnding = 'dim='+str(numVars)+'_deg='+str(degreeBound)+'_random'
    with open('data_'+fileEnding,'w') as f:
        f.write('bezout,LPRC,MV,gpSize\n')
        f.write('\n'.join(data))

#import cProfile
#cProfile.run('asdf','runStats')
