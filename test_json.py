from pathlib import Path
from pprint import pprint
import os
import sys
import json

from acados_template import AcadosOcp

def test_json(file_name: str | Path):
    ocp = AcadosOcp.from_json(file_name)
    
    ocp.json_file = str(ocp.json_file).replace('.json', '_out.json')
    ocp.make_consistent()
    ocp_dict = ocp.to_dict()
    pprint(ocp_dict)
    
    ocp.dump_to_json()


if __name__ == "__main__":
    sys_file = Path(__file__).parent / 'acados_ocp.json'
    test_json(sys_file)