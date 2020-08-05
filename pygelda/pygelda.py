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

import logging

import numpy as np

import _gelda

logger = logging.getLogger(__name__)


class Gelda(object):
    def __init__(self, edif, adif, fdif, neq, ndif):
        self.__edif = edif
        self.__adif = adif
        self.__fdif = fdif
        self.__neq = neq
        self.__ndif = ndif

    def e(self, t, idif, e, ipar, rpar, n, lde):
        # import ipdb; ipdb.set_trace()
        if idif > self.__ndif:
            ierr = -2
            return ierr
        else:
            e[:n, :n] = self.__edif(t, idif)
            ierr = 0
            return ierr

    def a(self, t, idif, a, ipar, rpar, n, lda):
        if idif > self.__ndif:
            ierr = -2
            return ierr
        else:
            a[:n, :n] = self.__adif(t, idif)
            ierr = 0
            return ierr

    def f(self, t, idif, f, ipar, rpar, n):
        if idif > self.__ndif:
            ierr = -2
            return ierr
        else:
            f[:] = self.__fdif(t, idif)
            ierr = 0
            return ierr

    def solve(self, t, x0, rtol=1e-6, atol=1e-6):
        # if x0.ndim != 1:
        #     raise ValueError("x0 needs 1d vector")
        lent = t.size
        if lent == 0:
            return np.zeros((0, )), 0
        xout = np.zeros((self.__neq, lent))
        rtol = np.ones((self.__neq, )) * rtol
        atol = np.ones((self.__neq, )) * atol
        ipar = np.arange(1, 2)
        rpar = 1. * ipar
        info_array = np.zeros((20, ))
        info_array[2] = 0
        t0 = np.array([t[0]])
        x0_f = x0.copy()
        for i, tn in enumerate(t):
            tout = np.array([tn])
            ierr = 1
            logger.debug("t0: {}\ntout: {}\nx0: {}".format(t0, tout, x0))
            while ierr == 1:
                xprime, cval, iwarn, ierr = _gelda.dgelda(self.e,
                                                          self.a,
                                                          self.f,
                                                          t0,
                                                          tout,
                                                          x0_f,
                                                          ipar,
                                                          rpar,
                                                          rtol,
                                                          atol,
                                                          info=info_array)
                xout[:, i] = x0_f
        return xout, ierr
