import unittest
import sys
sys.path.append("../module_py")
import dfm2

class TestCad2D(unittest.TestCase):
    def test_cad(self):
    	cad = dfm2.Cad2D()
    	cad.add_polygon([-1,-1, +1,-1, +1,+1, -1,+1])
    	cad = dfm2.Cad2D(list_xy=[-1,-1, +1,-1, +1,+1, -1,+1])
    	msh = cad.mesh(0.02)


if __name__ == "__main__":
    unittest.main()