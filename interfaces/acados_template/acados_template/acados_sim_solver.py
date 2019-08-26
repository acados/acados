from ctypes import *
import numpy as np

class acados_sim_solver:
    def __init__(self, acados_sim, shared_lib):

        self.sim_struct = acados_sim
        model_name = self.sim_struct.model_name


        self.shared_lib = CDLL(shared_lib)
        getattr(self.shared_lib, f"{model_name}_acados_sim_create")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_opts").restype = c_void_p
        self.sim_opts = getattr(self.shared_lib, f"{model_name}_acados_get_sim_opts")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_dims").restype = c_void_p
        self.sim_dims = getattr(self.shared_lib, f"{model_name}_acados_get_sim_dims")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_config").restype = c_void_p
        self.sim_config = getattr(self.shared_lib, f"{model_name}_acados_get_sim_config")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_out").restype = c_void_p
        self.sim_out = getattr(self.shared_lib, f"{model_name}_acados_get_sim_out")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_in").restype = c_void_p
        self.sim_in = getattr(self.shared_lib, f"{model_name}_acados_get_sim_in")()

        nu = self.sim_struct.dims.nu
        nx = self.sim_struct.dims.nx
        self.gettable = {
            'x': nx,
            'xn': nx,
            'u': nu,
            'S_forw': nx*(nx+nu)
        }

        self.settable = ['S_adj', 'S_forw', 'T', 'x', 'u', 'xdot', 'z', 'Su', 'Sx']
        self.model_name = model_name

    def solve(self):
        status = getattr(self.shared_lib, f"{self.model_name}_acados_sim_solve")()
        return status

    def get(self, field_):

        field = field_
        field = field.encode('utf-8')

        #  ocp_sim_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        #  self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int
        #  dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, field)

        if field_ in self.gettable.keys():

            # allocate array
            dims = self.gettable[field_]
            out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))

            self.shared_lib.sim_out_get.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
            self.shared_lib.sim_out_get(self.sim_config, self.sim_dims, self.sim_out, field, out_data);

        else:
            raise Exception(f'acados_solver.set(): Unknown field {field}, available fiels are {",".join(self.gettable.keys())}')

        # out = cast((out), POINTER(c_double))

        return out

    def set(self, field_, value_):

        # cast value_ to avoid conversion issues
        if type(value_) == float:
            value_ = np.array([value_])

        value_ = value_.astype(float)
        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)


        field = field_
        field = field.encode('utf-8')

        #  self.shared_lib.ocp_sim_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        #  self.shared_lib.ocp_sim_dims_get_from_attr.restype = c_int

        #  dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, field)

        #  if value_.shape[0] != dims:
            #  raise Exception('acados_solver.set(): mismatching dimension for field "{}" with dimension {} (you have {})'.format(field_,dims, value_.shape[0]))

        if field_ in self.settable:
            self.shared_lib.sim_in_set.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
            self.shared_lib.sim_in_set(self.sim_config, self.sim_dims, self.sim_in, field, value_data_p);
        else:
            raise Exception(f'acados_solver.set(): Unknown field {field}, available fiels are {",".join(self.settable)}')

        return

    def __del__(self):
        getattr(self.shared_lib, f"{self.model_name}_acados_sim_free")()
