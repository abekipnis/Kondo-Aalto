from AaltoAtoms import Spec
import numpy as np

test_file = r"tests\Createc2_210811.113827.L0001.VERT"

class TestSpec():

    def test_init(self):
        global test_file
        S = Spec(test_file)
        assert(S.fname==test_file)
        assert(S.NPoints==512)
        assert(np.isclose(S.FBLogiset, 1500, 0.1))
        assert(np.isclose(S.biasVoltage, 100, 1))
        assert(len(S.current)==len(S.dIdV))
