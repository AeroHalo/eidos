# Contributing to EIDOS

We welcome contributions to the EIDOS framework! This document provides guidelines for contributing.

## Code of Conduct

We adhere to a code of conduct that fosters an open and welcoming environment. Please be respectful and constructive in all interactions.

## How to Contribute

### Reporting Issues

1. Check existing issues to avoid duplicates
2. Use issue templates when available
3. Provide minimal reproducible examples
4. Include system information (Python version, OS, etc.)

### Submitting Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass: `pytest`
6. Format code: `black eidos/`
7. Submit PR with clear description

### Code Style

- Follow PEP 8
- Use type hints where appropriate
- Add docstrings to all functions/classes
- Keep line length under 88 characters (Black default)

### Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=eidos

# Run specific test
pytest tests/test_core.py::test_entropy_field
```
### Documentation

Update docstrings for API changes

Add examples to docstrings

Update tutorials if needed

Build docs locally: cd docs && make html

## Development Setup
```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/EIDOS.git
cd EIDOS

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows

# Install in development mode
pip install -e .[dev]

# Set up pre-commit hooks (optional)
pre-commit install
```
## Release Process

Update version in setup.py and eidos/__init__.py

Update CHANGELOG.md

Create release PR

Tag release: git tag -a v1.0.0 -m "Release version 1.0.0"

Push tags: git push origin --tags

## Questions?
Feel free to open a discussion or contact the maintainers!