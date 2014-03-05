#/bin/env python
import numpy as N

def gcd(a,b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:      
        a, b = b, a % b
    return a

class Frac(object):
    def __init__(self, num, denom):
        self.num = num
        self.denom = denom
        self.reduce()
        return
    
    def reduce(self):
        g = gcd(self.num, self.denom)
        self.num /= g
        self.denom /= g
        return
    
    def __add__(self, other):
        if type(other) is Frac:
            return Frac(self.num*other.denom + other.num*self.denom,
                       self.denom*other.denom)
        elif type(other) is int:
            return Frac(self.num + self.denom*other, self.denom)
        elif type(other) is float:
            return other + float(self.num) / float(self.denom)
        else:
            return NotImplemented
    def __radd__(self, other):
        return self+other
    def __sub__(self, other):
        return self + -other
    def __rsub__(self, other):
        return other +-self
    
    def __mul__(self, other):
        if type(other) is Frac:
            return Frac(self.num*other.num,
                       self.denom*other.denom)
        elif type(other) is int:
            return Frac(self.num*other, self.denom)
        elif type(other) is float:
            return (self.num* other)/ self.denom
        else:
            return NotImplemented
    def __rmul__(self, other):
        return self*other
    
    def __neg__(self):
        return Frac(-self.num, self.denom)
    
    def __str__(self):
        if self.num:
            return "%d./%d." % (self.num, self.denom)
        else:
            return "0"

    def __repr__(self):
        return "Frac(%d, %d)" % (self.num, self.denom)

    def __float__(self):
        return float(self.num)/float(self.denom)
    
    
ndim = 3
nvel = 15

cs2 = Frac(1,3)

# if you reorder these, must similarly reorder chi1
velocityList = [
    #velocity, weight, chi1
    ([ 0, 0, 0], Frac(2, 9), -2),
    ([ 1, 0, 0], Frac(1, 9),  1),
    ([-1, 0, 0], Frac(1, 9),  1),
    ([ 0, 1, 0], Frac(1, 9),  1),
    ([ 0,-1, 0], Frac(1, 9),  1),
    ([ 0, 0, 1], Frac(1, 9),  1),
    ([ 0, 0,-1], Frac(1, 9),  1),
    ([ 1, 1, 1], Frac(1,72), -2),
    ([ 1, 1,-1], Frac(1,72), -2),
    ([ 1,-1, 1], Frac(1,72), -2),
    ([ 1,-1,-1], Frac(1,72), -2),
    ([-1, 1, 1], Frac(1,72), -2),
    ([-1, 1,-1], Frac(1,72), -2),
    ([-1,-1, 1], Frac(1,72), -2),
    ([-1,-1,-1], Frac(1,72), -2),
    ]

xi = N.array([vl[0] for vl in velocityList])
w = N.array([vl[1] for vl in velocityList])

complement = N.zeros(nvel, dtype=N.int)
for i in range(nvel):
    for j in range(i, nvel):
        if N.alltrue(xi[i] == -xi[j]):
            complement[i] = j
            complement[j] = i
            break

Q = N.zeros((nvel, ndim, ndim), dtype=N.object)
for i in range(nvel):
    Q[i] = N.outer(xi[i], xi[i]) - N.identity(ndim, dtype=N.int) * cs2

chi1 = N.array([vl[2] for vl in velocityList])
jchi1 = chi1[:,N.newaxis]*xi
chi2 = N.multiply.reduce(xi, axis=-1)



dims = [('X', 0),
        ('Y', 1),
        ('Z', 2)]
       
modeList = [
    # name,     components,      normalizer
    ('rho',    N.ones(nvel,dtype=N.int), 1),
    ('momX',   xi[:,0],                  3),
    ('momY',   xi[:,1],                  3),
    ('momZ',   xi[:,2],                  3),
    ('SXX',    Q[:,0,0],         Frac(9,2)),
    ('SXY',    Q[:,0,1],                 9),
    ('SXZ',    Q[:,0,2],                 9),
    ('SYY',    Q[:,1,1],         Frac(9,2)),
    ('SYZ',    Q[:,1,2],                 9),
    ('SZZ',    Q[:,2,2],         Frac(9,2)),
    ('chi1',   chi1,             Frac(1,2)),
    ('jchi1X', jchi1[:,0],       Frac(3,2)),
    ('jchi1Y', jchi1[:,1],       Frac(3,2)),
    ('jchi1Z', jchi1[:,2],       Frac(3,2)),
    ('chi2',   chi2,                     9),
    ]

mm = N.zeros((nvel, nvel), dtype=N.object)
for i in range(nvel):
    mm[i] =modeList[i][1]
    
norms = N.array([ml[2] for ml in modeList])

mmi = (mm*w).transpose()*norms

for i in range(nvel):
    for j in range(nvel):
        if mmi[i,j].denom == 1:
            mmi[i,j] = mmi[i,j].num


def generateC():
    # generate some C now
    initCode = """
/* speed of sound */
lat->cs2 = %s;

/* weights of velocities */
""" % cs2

    for i in range (nvel):
        initCode += "lat->w[%d] = %s;\n" % (i, w[i])

    initCode += """
/* velocity vectors */
"""
    for i in range(nvel):
        for j in range(ndim):
            initCode += "lat->xi[%d][%d] = %s; " %(i,j, xi[i,j])
        initCode += "\n"

    initCode += """
/* lookup for complementary velocities */
"""
    for i in range(nvel):
        initCode += "lat->complement[%d] = %d;\n" % (i, complement[i])

    initCode += """
/* set up the Q_p_ab's */
"""
    for i in range(nvel):
        for a in range(ndim):
            for b in range(ndim):
                initCode += "lat->Q[%d][%d][%d] = %s; """ % (i,a,b, Q[i,a,b])
            initCode += "\n"
        initCode += "\n"

    initCode += """
/* Normalizers for modes */
"""
    for i in range(nvel):
        initCode += "lat->norms[%d] = %s;\n" % (i,norms[i])

    initCode += """
/* Matrix of eigenvectors. */
"""
    for i in range(nvel):
        for j in range(nvel):
            initCode += "lat->mm[%d][%d] = %s; " % (i,j, mm[i,j])
        initCode += "\n"

    initCode += """
/* And the inverse matrix... */
"""
    for i in range(nvel):
        for j in range(nvel):
            initCode += "lat->mmi[%d][%d] = %s; " % (i,j, mmi[i,j])
        initCode += "\n"

    # set include guards around the initCode
    initCode = "#ifdef DQ_init\n" + initCode + "#endif"

    # set up macro code
    macroCode = """
/* dimensions */
"""
    for d in dims:
        macroCode += '#define DQ_%s %d\n' % (d[0], d[1])

    macroCode += """
/* mode names */
"""
    for i in range(nvel):
        macroCode += '#define DQ_%s %d\n' % (modeList[i][0], i)

    def getMomMode(dInd):
        momName = 'mom'+[d for d in dims if d[1]==dInd][0][0]
        for i in range(nvel):
            if modeList[i][0] == momName:
                return i

    macroCode += "#define DQ_mom(dim) %d+dim\n" % getMomMode(0)

    macroCode = "#ifdef DQ_access_macros\n" + macroCode + "#endif\n\n"

    print "/* Generated by eigenvectors.py - edit that and rerun if needed */\n"
    print macroCode
    print initCode
    return


def generatePy():
    # generate some Python now
    code = ''
    
    code += '# Weights of each f_i\n'
    code += 'lat.w = N.array([' + \
            ', '.join([str(el) for el in w])+ \
            '], dtype=N.float)\n\n'
    
    code += '# Velocity vectors\n'
    code += 'lat.xi = N.array([' + \
            ',\n                  '.join(['['+', '.join([str(xiia) for xiia in xii]) +']'
                        for xii in xi]) + '], dtype=N.float)\n\n'

    code += '# complement of each velocity (xi[complement[i]] == -xi[i])\n'
    code += 'lat.complement = N.array([' + \
            ', '.join([str(ci) for ci in complement]) + \
            '], dtype=N.int)\n\n'
    
    code += '# Kinetic projectors\n'
    code += 'lat.Q = N.array([' + ',\n                 '.join(['['+ ',\n                  '.join(['[' + ', '.join([str(Qiab) for Qiab in Qia]) + ']' for Qia in Qi]) +']' for Qi in Q]) + '], dtype=N.float)\n\n'
    
    code += '# Normalizers for mm\n'
    code += 'lat.norms = N.array([' + \
            ', '.join([str(ni) for ni in norms]) + \
            '], dtype=N.float)\n\n'
    
    code += '# Mode matrix that converts from the velocity basis to the\n# moment basis (density, momentum, stress, ghosts, etc.;\n# ordering as specified in libd3q15/eigh) (i.e. m_i = \n# mm_ij f_j)\n'
    code += 'lat.mm = N.array([' + \
            ',\n                  '.join(['[' + ', '.join([str(mmab) for mmab in mma]) + ']' for mma in mm]) + \
    '], dtype=N.float)\n\n'
    
    code += '# Inverse of mm\n'
    code += 'lat.mmi = (lat.mm*lat.w).transpose()*lat.norms\n\n'

    code += '# the mode ordering\n'
    code += 'lat.modes = [' + \
        ', \n             '.join(['"%s"' % x[0] for x in modeList]) + \
        ']'
    
    # indent by 4 spaces
    code = """
import numpy as N
def init(lat):
    """ + code.replace('\n', '\n    ')
    
    print "# Generated by eigenvectors.py - edit that and rerun if needed \n"
    print code

if __name__ == "__main__":
    import sys
    
    if sys.argv[1].upper() == 'C':
        generateC()
    elif sys.argv[1].upper() == 'PY':
        generatePy()
    
