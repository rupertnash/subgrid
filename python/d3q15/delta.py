"""Simple python implementation of Peskin's delta function.

"""
from numpy import abs, sqrt

def delta(x):
    """Discretised Dirac delta function.

    Returns an approximation to the delta function, as explained in
    Peskin Acta Numerica (2002) pp. 479-517, chapter 6.
    """
    abs_x = abs(x)
    if abs_x <= 1.0:
        return (3 - 2*abs_x + sqrt(1 + 4*abs_x - 4*abs_x*abs_x))/8.
    elif abs_x <= 2.0:
        return (5 - 2*abs_x-sqrt(-7 + 12*abs_x - 4*abs_x*abs_x))/8.
    else:
        return 0.0
