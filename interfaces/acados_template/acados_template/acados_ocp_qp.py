import json
import numpy as np
from .utils import cast_to_2d_nparray, check_if_nparray_and_flatten, is_empty

class AcadosOcpQpDims:
    """
    Class containing the dimensions of an OCP-structured QP problem
    """
    def __init__(self, N: int):
        self.N = N
        self.nx = np.zeros((N + 1,), dtype=int)
        self.nu = np.zeros((N + 1,), dtype=int)
        self.nbx = np.zeros((N + 1,), dtype=int)
        self.nbu = np.zeros((N + 1,), dtype=int)
        self.nb = np.zeros((N + 1,), dtype=int)
        self.ng = np.zeros((N + 1,), dtype=int)
        self.ns = np.zeros((N + 1,), dtype=int)


class AcadosOcpQp:
    """
    Class containing the description of an OCP-structured QP problem
    """
    def __init__(self, N: int):
        self.__N = N
        # dynamics
        self.__A_list = [None] * N
        self.__B_list = [None] * N
        self.__b_list = [None] * N
        # cost
        self.__Q_list = [None] * (N + 1)
        self.__R_list = [None] * (N + 1)
        self.__S_list = [None] * (N + 1)
        self.__q_list = [None] * (N + 1)
        self.__r_list = [None] * (N + 1)
        self.__zl_list = [None] * (N + 1)
        self.__zu_list = [None] * (N + 1)
        self.__Zl_list = [None] * (N + 1)
        self.__Zu_list = [None] * (N + 1)
        ## constraints
        # bounds
        self.__idxb_list = [None] * (N + 1)
        self.__lbu_list = [None] * (N + 1)
        self.__ubu_list = [None] * (N + 1)
        self.__lbx_list = [None] * (N + 1)
        self.__ubx_list = [None] * (N + 1)
        # general linear constraints
        self.__C_list = [None] * (N + 1)
        self.__D_list = [None] * (N + 1)
        self.__lg_list = [None] * (N + 1)
        self.__ug_list = [None] * (N + 1)
        # slacks
        self.__idxs_rev_list = [None] * (N + 1)
        self.__lls_list = [None] * (N + 1)
        self.__lus_list = [None] * (N + 1)
        # masks
        self.__lbu_mask_list = [None] * (N + 1)
        self.__ubu_mask_list = [None] * (N + 1)
        self.__lbx_mask_list = [None] * (N + 1)
        self.__ubx_mask_list = [None] * (N + 1)
        self.__lg_mask_list = [None] * (N + 1)
        self.__ug_mask_list = [None] * (N + 1)
        self.__lls_mask_list = [None] * (N + 1)
        self.__lus_mask_list = [None] * (N + 1)

        self.__dims = AcadosOcpQpDims(N)

        # meta
        self.dynamics_fields = ['A', 'B', 'b']
        self.cost_fields = ['Q', 'R', 'S', 'q', 'r', 'zl', 'zu', 'Zl', 'Zu']
        self.constraint_fields = ['idxb', 'lbu', 'ubu', 'lbx', 'ubx',
                                  'C', 'D', 'lg', 'ug',
                                  'idxs_rev', 'lls', 'lus',
                                  'lbu_mask', 'ubu_mask', 'lbx_mask', 'ubx_mask',
                                  'lg_mask', 'ug_mask', 'lls_mask', 'lus_mask']
        self.all_fields = self.dynamics_fields + self.cost_fields + self.constraint_fields

        self.vector_fields = ['b', 'q', 'r', 'zl', 'zu', 'Zl', 'Zu',
                              'idxb', 'lbu', 'ubu', 'lbx', 'ubx',
                              'lg', 'ug',
                              'idxs_rev', 'lls', 'lus',
                              'lbu_mask', 'ubu_mask', 'lbx_mask', 'ubx_mask',
                              'lg_mask', 'ug_mask', 'lls_mask', 'lus_mask']
        self.matrix_fields = ['A', 'B', 'Q', 'R', 'S', 'C', 'D']

    @property
    def N(self) -> int:
        return self.__N

    @property
    def dims(self) -> AcadosOcpQpDims:
        return self.__dims

    # Dynamics properties
    @property
    def A(self) -> list:
        return self.__A_list

    @property
    def B(self) -> list:
        return self.__B_list

    @property
    def b(self) -> list:
        return self.__b_list

    # Cost properties
    @property
    def Q(self) -> list:
        return self.__Q_list

    @property
    def R(self) -> list:
        return self.__R_list

    @property
    def S(self) -> list:
        return self.__S_list

    @property
    def q(self) -> list:
        return self.__q_list

    @property
    def r(self) -> list:
        return self.__r_list

    @property
    def zl(self) -> list:
        return self.__zl_list

    @property
    def zu(self) -> list:
        return self.__zu_list

    @property
    def Zl(self) -> list:
        return self.__Zl_list

    @property
    def Zu(self) -> list:
        return self.__Zu_list

    # Constraint properties
    @property
    def idxb(self) -> list:
        return self.__idxb_list

    @property
    def lbu(self) -> list:
        return self.__lbu_list

    @property
    def ubu(self) -> list:
        return self.__ubu_list

    @property
    def lbx(self) -> list:
        return self.__lbx_list

    @property
    def ubx(self) -> list:
        return self.__ubx_list

    @property
    def C(self) -> list:
        return self.__C_list

    @property
    def D(self) -> list:
        return self.__D_list

    @property
    def lg(self) -> list:
        return self.__lg_list

    @property
    def ug(self) -> list:
        return self.__ug_list

    @property
    def idxs_rev(self) -> list:
        return self.__idxs_rev_list

    @property
    def lls(self) -> list:
        return self.__lls_list

    @property
    def lus(self) -> list:
        return self.__lus_list

    # Mask properties
    @property
    def lbu_mask(self) -> list:
        return self.__lbu_mask_list

    @property
    def ubu_mask(self) -> list:
        return self.__ubu_mask_list

    @property
    def lbx_mask(self) -> list:
        return self.__lbx_mask_list

    @property
    def ubx_mask(self) -> list:
        return self.__ubx_mask_list

    @property
    def lg_mask(self) -> list:
        return self.__lg_mask_list

    @property
    def ug_mask(self) -> list:
        return self.__ug_mask_list

    @property
    def lls_mask(self) -> list:
        return self.__lls_mask_list

    @property
    def lus_mask(self) -> list:
        return self.__lus_mask_list

    def set(self, field_name: str, stage: int, value: np.ndarray):
        if stage < 0 or stage > self.N:
            raise ValueError(f"Stage {stage} is out of bounds for N={self.N}.")
        elif field_name in self.dynamics_fields and stage == self.N:
            raise ValueError(f"Dynamics fields cannot be set at terminal stage N={self.N}.")

        if field_name in self.all_fields:
            field_list = getattr(self, f'_{self.__class__.__name__}__{field_name}_list')
            if field_name in self.vector_fields:
                field_list[stage] = check_if_nparray_and_flatten(value, name=f"{field_name} at stage {stage}")
            elif field_name in self.matrix_fields:
                field_list[stage] = cast_to_2d_nparray(value, name=f"{field_name} at stage {stage}")
        else:
            raise ValueError(f"Field name {field_name} is not recognized.")

    def make_consistent(self, assert_dims: bool = True):
        # detect dims
        nx_next = None

        for i in range(self.N+1):
            # cost
            nx = self.Q[i].shape[0]
            nu = self.R[i].shape[0]
            self.__dims.nx[i] = nx
            self.__dims.nu[i] = nu

            if assert_dims:
                assert self.q[i].shape == (nx,), f"Inconsistent dimensions in q vector at stage {i}."
                assert self.r[i].shape == (nu,), f"Inconsistent dimensions in r vector at stage {i}."
                if nu > 0:
                    assert self.S[i].shape == (nu, nx), f"Inconsistent dimensions in S matrix at stage {i}."

            # dynamics
            if assert_dims and i < self.N:
                if nx_next is not None:
                    assert nx == nx_next, f"Inconsistent dimensions between consecutive A matrices at stage {i}."

                nx_next = self.A[i].shape[0]

                assert self.B[i].shape[0] == nx_next, "Inconsistent dimensions between A and B matrices."
                assert self.b[i].shape[0] == nx_next, "Inconsistent dimensions between A and b matrices."

            # constraints
            self.__dims.nbx[i] = len(self.lbx[i])
            self.__dims.nbu[i] = len(self.lbu[i])
            self.__dims.nb[i] = self.__dims.nbx[i] + self.__dims.nbu[i]
            self.__dims.ng[i] = len(self.lg[i])
            self.__dims.ns[i] = len(self.lls[i])

            if assert_dims:
                assert len(self.idxb[i]) == self.__dims.nb[i], f"Inconsistent number of bound constraint indices at stage {i}."
                assert self.C[i].shape[0] == self.__dims.ng[i], f"Inconsistent number of general constraints at stage {i}."
                assert self.D[i].shape[0] == self.__dims.ng[i], f"Inconsistent number of general constraints at stage {i}."
                if self.__dims.ng[i] > 0:
                    assert self.C[i].shape[1] == self.__dims.nx[i], f"Inconsistent number of states in general constraints at stage {i}."
                    assert self.D[i].shape[1] == self.__dims.nu[i], f"Inconsistent number of inputs in general constraints at stage {i}."

                assert len(self.idxs_rev[i]) == self.__dims.nb[i] + self.__dims.ng[i], f"Inconsistent number of slack variable indices at stage {i}."

                # slack cost
                assert self.zl[i].shape == (self.__dims.ns[i],), f"Inconsistent dimensions in zl vector at stage {i}."
                assert self.zu[i].shape == (self.__dims.ns[i],), f"Inconsistent dimensions in zu vector at stage {i}."
                assert self.Zl[i].shape == (self.__dims.ns[i],), f"Inconsistent dimensions in Zl vector at stage {i}."
                assert self.Zu[i].shape == (self.__dims.ns[i],), f"Inconsistent dimensions in Zu vector at stage {i}."

                # masks
                assert self.lbu_mask[i].shape == (self.__dims.nbu[i],), f"Inconsistent dimensions in lbu_mask at stage {i}."
                assert self.ubu_mask[i].shape == (self.__dims.nbu[i],), f"Inconsistent dimensions in ubu_mask at stage {i}."
                assert self.lbx_mask[i].shape == (self.__dims.nbx[i],), f"Inconsistent dimensions in lbx_mask at stage {i}."
                assert self.ubx_mask[i].shape == (self.__dims.nbx[i],), f"Inconsistent dimensions in ubx_mask at stage {i}."
                assert self.lg_mask[i].shape == (self.__dims.ng[i],), f"Inconsistent dimensions in lg_mask at stage {i}."
                assert self.ug_mask[i].shape == (self.__dims.ng[i],), f"Inconsistent dimensions in ug_mask at stage {i}."
                assert self.lls_mask[i].shape == (self.__dims.ns[i],), f"Inconsistent dimensions in lls_mask at stage {i}."
                assert self.lus_mask[i].shape == (self.__dims.ns[i],), f"Inconsistent dimensions in lus_mask at stage {i}."


    @classmethod
    def from_dict(cls, qp_dict) -> 'AcadosOcpQp':
        # Determine N from the keys
        Q_keys = [k for k in qp_dict.keys() if k.startswith('Q_')]
        N = len(Q_keys) - 1
        lN = len(str(N+1))

        # Create instance
        qp = cls(N)

        # Parse dynamics matrices (stages 0 to N-1)
        for field in qp.dynamics_fields:
            for i in range(N):
                key = f'{field}_{i:0{lN}d}'
                if key in qp_dict:
                    qp.set(field, i, qp_dict[key])
                else:
                    qp.set(field, i, np.zeros((0,0)))

        # Parse cost & constraint data (stages 0 to N)
        for field in qp.cost_fields + qp.constraint_fields:
            for i in range(N + 1):
                key = f'{field}_{i:0{lN}d}'
                if key in qp_dict:
                    qp.set(field, i, qp_dict[key])
                else:
                    qp.set(field, i, np.zeros((0,0)))

        qp.make_consistent()
        return qp

    @classmethod
    def from_json(cls, json_file: str) -> 'AcadosOcpQp':
        with open(json_file, 'r') as f:
            qp_dict = json.load(f)

        # Convert lists to numpy arrays
        for key, value in qp_dict.items():
            if isinstance(value, list):
                if any(field in key for field in ['idxb', 'idxs_rev']):
                    qp_dict[key] = np.array(value, dtype=int)
                else:
                    qp_dict[key] = np.array(value)

        return cls.from_dict(qp_dict)
