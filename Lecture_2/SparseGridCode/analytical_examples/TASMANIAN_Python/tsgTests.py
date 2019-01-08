#!/usr/bin/python

import unittest
import TasmanianSG
import math
import numpy

grid = TasmanianSG.TasmanianSparseGrid()

class TestSomething(unittest.TestCase):
        
    def testGlobal(self):
        pass

    def test_example1(self):
        dim = 2
        for lvl, tol in zip([6, 7], [1e-11, 1e-14]):
            grid.makeGlobalGrid(dim, 0, lvl, "level", "clenshaw-curtis")
            points = grid.getPoints()
            weights = grid.getQuadratureWeights()
            integral = sum([weights[i] * math.exp(- points[i][0] * points[i][0]) * math.cos(points[i][1]) for i in range(points.shape[0])])
            err = math.fabs(integral - 2.513723354063905e+00)
            self.assertLessEqual(err, tol)

    def test_example2(self):
        dim = 2
        for prec, tol in zip([20, 40], [1e-2, 1e-10]):
            grid.makeGlobalGrid(dim, 0, prec, "qptotal", "gauss-patterson")
            grid.setDomainTransform(numpy.array([[-5.0, 5.0], [-2.0, 3.0]]))
            points = grid.getPoints()
            weights = grid.getQuadratureWeights()
            integral = sum([weights[i] * math.exp(- points[i][0] * points[i][0]) * math.cos(points[i][1]) for i in range(points.shape[0])])
            err = math.fabs(integral - 1.861816427518323e+00)
            self.assertLessEqual(err, tol)

    def test_raise_errors(self):
        return # comment this line to see what happen
        dim = 2
        lvl = 6
        grid.makeGlobalGrid(dim, 0, lvl, "level", "notregistered" )
        self.assertRaises(grid.makeGlobalGrid(dim, 0, lvl, "invalid", "clenshaw-curtis"), RuntimeError)

if __name__ == '__main__':
    unittest.main()

