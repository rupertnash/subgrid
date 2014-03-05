#!/usr/bin/env python
import sys
import os.path
from distutils.core import setup, Extension
import distutils.sysconfig
import numpy

numpyIncludeDir = numpy.get_include()

#compile_args = []
compile_args = ['-m32']
#link_args = []
link_args = ['-m32']
if os.path.exists('debug'):
    compile_args += ['-O0', '-U NDEBUG']
    
libSrcDir = '../libd3q15/'

cydelta_c = os.path.join('d3q15', 'CyDelta.c')
if not os.path.exists(cydelta_c):
    # try to generate it with Cython
    np_base = os.path.split(numpy.__file__)[0]
    cy_np_dir = os.path.join(np_base, 'doc', 'cython')
    cydelta_pyx = os.path.join('d3q15', 'CyDelta.pyx')
    
    cmd = 'cython -I%s -o %s %s' % (cy_np_dir,
                                    cydelta_c,
                                    cydelta_pyx)
    if os.system(cmd):
        raise EnvironmentError("Problem generating CyDelta.c with command: %s" % cmd)
    pass


ext_modules = [Extension('d3q15._d3q15', ['d3q15/d3q15.i',
                                          libSrcDir+'d3q15.c',
                                          libSrcDir+'noise.c',
                                          libSrcDir+'checks.c',
                                          libSrcDir+'bc_pbc.c',
                                          libSrcDir+'bc_noslip.c',
                                          libSrcDir+'bc_freeslip.c',
                                          libSrcDir+'bc_wall.c',
                                          libSrcDir+'force_none.c',
                                          libSrcDir+'force_const.c'],
                         include_dirs=[libSrcDir, numpyIncludeDir],
                         extra_compile_args=['-std=c99'] + compile_args,
                         extra_link_args=link_args),
               
               Extension('d3q15._CDelta', ['d3q15/CDelta.i',
                                           'd3q15/CDelta.c'],
                         include_dirs=[numpyIncludeDir,],
                         extra_compile_args=compile_args,
                         extra_link_args=link_args),
               
               Extension('d3q15.CyDelta', ['d3q15/CyDelta.c'],
                         include_dirs=[numpyIncludeDir],
                         extra_compile_args=compile_args,
                         extra_link_args=link_args),
               
               Extension('d3q15._utils', ['d3q15/utils.i',
                                          'd3q15/utils.c',
                                          'd3q15/CDelta.c'],
                         include_dirs=[libSrcDir, numpyIncludeDir],
                         extra_compile_args=compile_args,
                         extra_link_args=link_args),
               
               Extension('d3q15.fallers.__helpers',
                         ['d3q15/fallers/helpers.i',
                          'd3q15/fallers/_helpers.c',
                          libSrcDir+'noise.c',
                          'd3q15/XArray.c',
                          'd3q15/CDelta.c',
                          'd3q15/utils.c'],
                         include_dirs=[libSrcDir, numpyIncludeDir],
                         extra_compile_args=compile_args,
                         extra_link_args=link_args),
               ]


# run the script to generate eigenvectors.h
cmd = '%s %s C > %s' % (sys.executable,
                        os.path.join(os.path.dirname(__file__),
                                     libSrcDir, 'eigenvectors.py'),
                        os.path.join(os.path.dirname(__file__),
                                     libSrcDir, 'eigenvectors.h'))
if os.system(cmd):
    raise EnvironmentError("Problem generating eigenvectors.h with command: %s" % cmd)

# run the script to generate eigenvectors.py
cmd = '%s %s py > %s' % (sys.executable,
                      os.path.join(os.path.dirname(__file__), libSrcDir, 'eigenvectors.py'),
                      os.path.join(os.path.dirname(__file__), 'd3q15', 'eigenvectors.py'))
if os.system(cmd):
    raise EnvironmentError("Problem generating eigenvectors.py with command: %s" % cmd)

setup(name='d3q15',
      version='1.0',
      description='3 dimensional, 15 velocity Lattice Boltzmann code.',
      author='Rupert Nash',
      author_email='rupert.nash@ed.ac.uk',
      packages=['d3q15', 'd3q15.fallers'],
#      ext_package='d3q15',
      ext_modules=ext_modules
      )



setup(name='dqTools',
      version='1.0',
      description='Run time and analysis tools for use with d3q15',
      author='Rupert Nash',
      author_email='rupert.nash@ed.ac.uk',
      packages=['dqTools']
      )
