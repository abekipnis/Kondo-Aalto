from AaltoAtoms import CircCorralData

class TestCircCorralData:
    test_file = r"data\03-24\1\A220323.115455.dat"

    def test_init():
        C = CircCorralData(test_file, "label")
        assert(C.file==test_file)
        assert(C.label=="label")
        assert(C.imshape==(256, 256))
        assert(C.centroids==None)
