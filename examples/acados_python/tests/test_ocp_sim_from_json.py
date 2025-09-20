from pathlib import Path
import sys
import json
import os
import numpy as np
import scipy.linalg

sys.path.insert(0, '../pendulum_on_cart/common')
from pendulum_model import export_pendulum_ode_model

from acados_template import AcadosOcp, AcadosSim


def compare_dicts(original_dct: dict, copy_dct: dict):
    faulty_keys = set()
    keys1 = set(original_dct.keys())
    keys2 = set(copy_dct.keys())
    for key in keys1 - keys2:
        faulty_keys.add(key)
        print(f"Key {key} not in copied json file")
    for key in keys2 - keys1:
        faulty_keys.add(key)
        print(f"Key {key} not in original json file")

    for key in keys1 & keys2:
        if key == 'json_file':
            continue
        value1 = original_dct.get(key, None)
        value2 = copy_dct.get(key, None)
        if isinstance(value1, dict) and isinstance(value2, dict):
            faulty_keys.update(compare_dicts(value1, value2))
        elif isinstance(value1, list) and isinstance(value2, list):
            if len(value1) != len(value2):
                faulty_keys.add(key)
                print(f"List for key {key} differs in length between json files")
            else:
                for v1, v2 in zip(value1, value2):
                    if isinstance(v1, dict) and isinstance(v2, dict):
                        faulty_keys.update(compare_dicts(v1, v2))
                    elif v1 != v2:
                        faulty_keys.add(key)
                        print(f"Value for key {key} differs between json files: {v1} != {v2}")
        else:
            if original_dct[key] != copy_dct[key]:
                faulty_keys.add(key)
                print(f"Value for key {key} differs between json files: {original_dct[key]} != {copy_dct[key]}")
        
    return faulty_keys


def compare_jsons(json_file_1: str | Path, json_file_2: str | Path):
    with open(json_file_1, 'r') as f:
        json_1 = json.load(f)
    with open(json_file_2, 'r') as f:
        json_2 = json.load(f)

    faulty_keys = compare_dicts(json_1, json_2)
    if len(faulty_keys) > 0:
        raise Exception(f"Differences found between the two json files: {faulty_keys}")


def create_ocp_json(file_name: str | Path):
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N = 20

    # set dimensions
    ocp.solver_options.N_horizon = N

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))
    
    # bound on u
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # initial state
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.json_file = str(file_name)
    ocp.code_export_directory = os.path.abspath(ocp.code_export_directory)
    ocp.make_consistent()
    ocp.generate_external_functions()
    ocp.dump_to_json()
    
    
def create_sim_json(file_name: str | Path):
    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()

    # set simulation time
    sim.solver_options.T = 0.1
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    
    sim.make_consistent()
    sim.dump_to_json(file_name)


def test_ocp_json(file_name: str | Path):
    print(f"test_ocp_json with file {file_name}")
    create_ocp_json(file_name)
    ocp = AcadosOcp.from_json(file_name)
    
    copied_json_file = str(file_name).replace('.json', '_copy.json')
    ocp.json_file = copied_json_file
    ocp.make_consistent()
    ocp.generate_external_functions()
    ocp.dump_to_json()
    compare_jsons(file_name, copied_json_file)
    
    
def test_sim_json(file_name: str | Path):
    print(f"test_sim_json with file {file_name}")
    create_sim_json(file_name)
    sim = AcadosSim.from_json(file_name)
    
    copied_json_file = str(file_name).replace('.json', '_copy.json')
    sim.make_consistent()
    sim.dump_to_json(copied_json_file)
    compare_jsons(file_name, copied_json_file)
    


if __name__ == "__main__":
    sys_file = Path().home() / 'acados_ocp.json'
    test_ocp_json(sys_file)
    
    sys_file = Path().home() / 'acados_sim.json'
    test_sim_json(sys_file)
    