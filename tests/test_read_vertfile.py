from AaltoAtoms import Spec

class TestSpec():
    import numpy as np
    test_file = r"tests\Createc2_210811.113827.L0001.VERT"

    def test_init(self):
        S = Spec(test_file)
        assert(S.fname==test_file)
        assert(S.NPoints==512)
        assert(np.isclose(S.FBLogiset, 1500, 0.1))
        assert(np.isclose(S.biasVoltage, 100, 1))
        assert(len(S.current)==len(S.dIdV))
