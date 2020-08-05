#!/usr/bin/env python
##
## Copyright (c) 2020
##
## @author: Daniel Bankmann
## @company: Technische Universität Berlin
##
##
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info
import glob
from os.path import join


lapack_opt = get_info('lapack_opt')
_base_name = 'pygelda'
wrapper = Extension('_gelda', [join(_base_name, "pygelda.pyf")] + glob.glob("GELDA/SRC/*.f"),
                    extra_f77_compile_args=["-ffixed-form", "-llapack", "-lblas"],
                    extra_f90_compile_args=["-ffixed-form", "-llapack", "-lblas"],
                    f2py_options=[],
                    **lapack_opt)
setup(
    name=_base_name,
    version='0.1.0',
    description='Wrapper for the FORTRAN code GELDA. Solves linear time-varying differential-algebraic equations',
    author='Daniel Bankmann',
    author_email='bankmann@math.tu-berlin.de',
    license='GPL',
    platforms='Linux',
    packages=['pygelda'],
    ext_modules = [wrapper],
)
