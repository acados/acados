from pathlib import Path
from pprint import pprint
import os
import sys
import json

from acados_template import AcadosOcp, AcadosSim

def test_ocp_json(file_name: str | Path):
    ocp = AcadosOcp.from_json(file_name)
    
    ocp.json_file = str(ocp.json_file).replace('.json', '_out.json')
    ocp.make_consistent()
    ocp_dict = ocp.to_dict()
    pprint(ocp_dict)
    
    ocp.dump_to_json()
    
    
def test_sim_json(file_name: str | Path):
    sim = AcadosSim.from_json(file_name)
    
    json_file = str(file_name).replace('.json', '_out.json')
    sim.make_consistent()
    sim_dict = sim.to_dict()
    pprint(sim_dict)
    
    sim.dump_to_json(json_file)
    


if __name__ == "__main__":
    sys_file = Path().home() / 'acados_ocp.json'
    test_ocp_json(sys_file)
    
    sys_file = Path().home() / 'acados_sim.json'
    test_sim_json(sys_file)
    