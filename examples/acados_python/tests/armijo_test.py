from pathlib import Path


if __name__ == "__main__":
    import pytest

    raise SystemExit(pytest.main([str(Path(__file__).with_name("test_armijo.py"))]))
