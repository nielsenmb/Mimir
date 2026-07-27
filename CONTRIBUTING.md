# Contributing

Mimir is in its initial development phase. Bug reports, numerical test cases,
and focused pull requests are welcome.

## Development setup

```bash
python -m pip install -e ".[test,mast,docs]"
```

Run the complete local check suite before submitting a pull request:

```bash
ruff check .
pytest
python -m build
sphinx-build -W -b html docs docs/_build/html
```

## Conventions

- Use NumPy-style docstrings for all public classes and functions.
- Add regression tests for every bug fix.
- Keep network access out of the default test suite; mock MAST calls or mark
  explicit integration tests separately.
- Document any change to units, normalization, or numerical behaviour.
- Keep PBjam compatibility changes separate from deliberate numerical changes.
