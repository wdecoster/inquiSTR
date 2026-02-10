# Contributing to inquiSTR

Thank you for your interest in contributing to inquiSTR! This document provides guidelines and instructions for contributing to the project. If you are unsure about anything, please feel to open an issue before starting.

## Development Setup

```bash
# Clone the repository
git clone https://github.com/wdecoster/inquiSTR.git
cd inquiSTR

# Install development tools
make setup

# Install git hooks for automated code quality checks
make install-hooks
```

## Code Quality

This project uses automated code quality checks to maintain consistency and prevent CI failures:

- **Formatting**: `cargo fmt` automatically formats code according to Rust standards
- **Linting**: `cargo clippy` checks for common mistakes and suggests improvements
- **Testing**: `cargo test` runs the test suite

### Git Hooks (Recommended)

The project includes pre-commit and pre-push hooks that automatically run quality checks:

```bash
# Install hooks (run once after cloning)
make install-hooks

# The hooks will now automatically:
# - Format code on commit (pre-commit)
# - Run clippy and formatting checks before push (pre-push)
```

The most common cause of CI failures is code formatting and clippy warnings. The git hooks prevent these issues by:

1. **Pre-commit hook**: Automatically formats your code before each commit
2. **Pre-push hook**: Ensures code passes formatting and clippy checks before pushing

This saves time by catching issues locally instead of discovering them in CI.

### Manual Quality Checks

You can also run quality checks manually:

```bash
# Format code
make fmt

# Run clippy linter
make clippy

# Run all pre-push checks
make pre-push

# Run all CI checks (formatting, linting, tests)
make ci

# Test preset catalog URLs (requires network access)
cargo test test_preset_urls -- --ignored --nocapture
```

### Available Make Targets

```bash
make help          # Show all available targets
make fmt           # Format code with cargo fmt
make clippy        # Run clippy linter with warnings as errors
make test          # Run tests
make build         # Build in release mode
make musl          # Build static binary with musl
make ci            # Run all CI checks (fmt-check, clippy, test)
make pre-push      # Run pre-push checks (fmt, clippy)
make install-hooks # Install git hooks for automated checks
```

## How to Contribute

1. **Fork and clone** the repository
2. **Install hooks**: `make install-hooks`
3. **Make changes**, add tests, run tests and commit your code (hooks will auto-format)
4. **Push changes** (hooks will run quality checks)
5. **Create a pull request**

The automated hooks ensure your contributions meet the project's quality standards before they reach CI.

## Questions or Issues?

If you have questions or run into issues, please open an issue on GitHub. We're happy to help!
