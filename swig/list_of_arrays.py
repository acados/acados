from numpy import copyto

class list_of_arrays(list):
    def __init__(self):
        super().__init__()
    def __setitem__(self, k, v):
        copyto(self.__getitem__(k), v)
