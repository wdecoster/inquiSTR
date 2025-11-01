# Makefile for inquiSTR development

.PHONY: all build test clean fmt clippy audit docs install help install-hooks pre-push

# Default target
all: fmt clippy test build

# Build the project
build:
	cargo build --release

# Build a statically-linked Linux binary using MUSL
.PHONY: build-musl
build-musl:
	@echo "Building static MUSL binary (x86_64-unknown-linux-musl)"
	rustup target add x86_64-unknown-linux-musl >/dev/null 2>&1 || true
	@if command -v cross >/dev/null 2>&1; then \
		echo "Using cross for reproducible musl build"; \
		OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
		cross build --release --target x86_64-unknown-linux-musl; \
	else \
		echo "Using cargo. Ensure musl-gcc is available (sudo apt-get install musl-tools)"; \
		OPENSSL_STATIC=1 LIBZ_SYS_STATIC=1 BZIP2_STATIC=1 ZSTD_STATIC=1 LZMA_API_STATIC=1 CURL_STATIC=1 \
		cargo build --release --target x86_64-unknown-linux-musl; \
	fi
	@echo "Binary: target/x86_64-unknown-linux-musl/release/inquiSTR"

# Run tests  
test:
	cargo test

# Clean build artifacts
clean:
	cargo clean

# Format code
fmt:
	cargo fmt

# Check formatting
fmt-check:
	cargo fmt --check

# Run clippy
clippy:
	cargo clippy --all-targets --all-features -- -D warnings

# Security audit
audit:
	cargo audit

# Check for outdated dependencies
outdated:
	cargo outdated --root-deps-only

# Generate documentation
docs:
	cargo doc --no-deps --open

# Install locally
install:
	cargo install --path .

# Run all checks (CI simulation)
ci: fmt-check clippy test
	@echo "All CI checks passed!"

# Development setup
setup:
	rustup component add rustfmt clippy
	cargo install cargo-audit cargo-outdated

# Install git hooks for automated checks
install-hooks:
	@echo "Installing git hooks..."
	@if [ ! -d .git ]; then \
		echo "❌ Not a git repository (missing .git directory)."; \
		echo "   Run 'git init' or clone the repo with git to enable hooks."; \
		exit 1; \
	fi
	@mkdir -p .git/hooks
	@cp -f .githooks/pre-commit .git/hooks/pre-commit
	@cp -f .githooks/pre-push .git/hooks/pre-push
	@chmod +x .git/hooks/pre-commit .git/hooks/pre-push
	@echo "✅ Git hooks installed successfully!"
	@echo "💡 The hooks will now run automatically on commit and push"

# Run all pre-push checks manually
pre-push: fmt clippy
	@echo "🎉 All pre-push checks passed!"

# Benchmark (if benchmarks exist)
bench:
	cargo bench

# Check everything is ready for commit
pre-commit: fmt clippy test
	@echo "Ready for commit!"

# Show help
help:
	@echo "Available targets:"
	@echo "  all        - Format, lint, test, and build"
	@echo "  build      - Build the project in release mode" 
	@echo "  test       - Run tests"
	@echo "  clean      - Clean build artifacts"
	@echo "  fmt        - Format code"
	@echo "  fmt-check  - Check if code is formatted"
	@echo "  clippy     - Run clippy linter"
	@echo "  audit      - Run security audit"
	@echo "  outdated   - Check for outdated dependencies"
	@echo "  docs       - Generate and open documentation"
	@echo "  install    - Install inquiSTR locally"
	@echo "  ci         - Run all CI checks"
	@echo "  setup      - Install required tools"
	@echo "  install-hooks - Install git hooks for automated checks"
	@echo "  pre-push   - Run all pre-push checks manually"
	@echo "  bench      - Run benchmarks"
	@echo "  pre-commit - Check everything before committing"
	@echo "  help       - Show this help message"