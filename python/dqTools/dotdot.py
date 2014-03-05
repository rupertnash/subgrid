"""Add the parent directory of the current working dir to the python
module search path.
"""
import sys
import os.path

dotdot = os.path.abspath('..')
if dotdot not in sys.path:
    sys.path.append(dotdot)
    pass

