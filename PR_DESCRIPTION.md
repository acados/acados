# Python: add AcadosOcpQpSolver

This PR adds Python bindings for directly solving OCP QP (Optimal Control Problem Quadratic Programming) problems.

## Key Changes

- **New Python classes**:
  - `AcadosOcpQp`: Represents an OCP QP problem, can be loaded from JSON
  - `AcadosOcpQpSolver`: Solver for OCP QP problems

- **C interface updates**:
  - Enhanced `ocp_qp_interface.c/h` with additional getter/setter functions
  - Updated QP common functionality

- **Testing**:
  - Added test suite in `examples/acados_python/tests/qp_test/`
  - Tests verify QP solutions match SQP iteration solutions
  - Integrated into CI workflow

## Usage Example

```python
from acados_template import AcadosOcpQp, AcadosOcpQpSolver

# Load QP from JSON file
qp = AcadosOcpQp.from_json('last_qp.json')

# Create solver and solve
solver = AcadosOcpQpSolver(qp)
solver.solve()

# Get solution
sol = solver.get_iterate()
```

This enables users to extract and solve the QP subproblems from SQP iterations independently.
