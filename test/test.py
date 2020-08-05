from pygelda import Gelda
import numpy as np

def edif(t, idif):
    e = np.diag([1,1])
    return e

def adif(t, idif):
    a = np.diag([1,2])
    return a

def fdif(t, idif):
    f = np.array([1.,0.])
    return f
def test_d0():
    test = Gelda(edif,adif,fdif,2,0)
    t = np.array([1.,2.,3.])
    x0 = np.array([1.,1.])
    xout, ierr = test.solve(t, x0)
    xout_ref = np.array([[1., 1., 1.],
       [1., 1., 1.]])
    assert np.allclose(xout,xout_ref)


def edif(t, idif):
    e = (1-idif)*np.diag([1,0])
    return e

def adif(t, idif):
    a = (1-idif)*np.array([[1,1],[1,0]])
    return a

def fdif(t, idif):
    if idif==0:
        f =  np.array([t,t])
    else:
        f =  np.array([1.,1.])
    return f
def test_s1():
        test2 = Gelda(edif,adif,fdif,2,1)
        t = np.array([1.,2.,3.])
        x0 = np.array([1.,-1.])
        xout, ierr = test2.solve(t, x0)
        xout_ref = np.array([[-1., -2., -3.], [-1., -1., -1.]])
        assert np.allclose(xout,xout_ref)
