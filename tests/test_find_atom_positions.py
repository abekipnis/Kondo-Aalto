from AaltoAtoms import CircCorralData
test_file = r"data\03-24\1\A220323.115455.dat"

class TestCircCorralData:

    def test_init(self):
        global test_file
        C = CircCorralData(test_file, "label")
        assert(C.file==test_file)
        assert(C.label=="label")
        assert(C.imshape==(256, 256))
        assert(C.centroids==None)
