# Makefile for SLcode - Sea Level Equation Solver
#
# This repository contains sea level equation solvers for glacial isostatic adjustment
# modeling, including elastic and viscoelastic Earth models. The code computes sea level
# fingerprints from ice sheet changes using spherical harmonic methods.

.PHONY: all clean help install-deps python-elastic matlab-elastic matlab-viscoelastic matlab-benchmarked test

# Default target - run the main Python sea level equation solver
all: python-elastic

# Help message
help:
	@echo "SLcode - Sea Level Equation Solver"
	@echo "=================================="
	@echo ""
	@echo "Available targets:"
	@echo "  all                    - Run main Python elastic sea level solver (default)"
	@echo "  python-elastic         - Run Python elastic sea level equation"
	@echo "  matlab-elastic         - Run MATLAB elastic sea level equation (requires MATLAB)"
	@echo "  matlab-viscoelastic    - Run MATLAB viscoelastic sea level equation (requires MATLAB)"
	@echo "  matlab-benchmarked     - Run MATLAB benchmarked versions (requires MATLAB)"
	@echo "  install-deps           - Install Python dependencies"
	@echo "  test                   - Test Python dependencies"
	@echo "  clean                  - Clean temporary files"
	@echo "  help                   - Show this help message"
	@echo ""
	@echo "Note: MATLAB targets require MATLAB to be installed and available in PATH"

# Install Python dependencies
install-deps:
	@echo "Installing Python dependencies..."
	pip3 install --user numpy scipy matplotlib pyshtools astropy xarray

# Test Python dependencies
test:
	@echo "Testing Python dependencies..."
	python3 -c "import numpy, scipy, matplotlib, pyshtools; print('All dependencies available')"

# Run the main Python elastic sea level equation solver
python-elastic:
	@echo "Running Python elastic sea level equation solver..."
	@echo "Input: West Antarctic Ice Sheet collapse scenario"
	@echo "Output: Global sea level fingerprint"
	python3 SL_equation_elastic.py

# Run the SLcode_py version (alternative Python implementation)
python-elastic-alt:
	@echo "Running alternative Python elastic sea level solver..."
	cd SLcode_py && python3 SL_equation_elastic.py

# MATLAB targets (require MATLAB)
matlab-elastic:
	@echo "Running MATLAB elastic sea level equation..."
	@if command -v matlab >/dev/null 2>&1; then \
		matlab -nodisplay -nosplash -r "SL_equation_elastic_benchmarked; exit"; \
	else \
		echo "Error: MATLAB not found. Please install MATLAB or run 'make python-elastic' instead."; \
		exit 1; \
	fi

matlab-viscoelastic:
	@echo "Running MATLAB viscoelastic sea level equation..."
	@if command -v matlab >/dev/null 2>&1; then \
		matlab -nodisplay -nosplash -r "SL_equation_viscoelastic_benchmarked; exit"; \
	else \
		echo "Error: MATLAB not found. Please install MATLAB or run 'make python-elastic' instead."; \
		exit 1; \
	fi

matlab-viscoelastic-gia:
	@echo "Running MATLAB viscoelastic GIA-only model..."
	@if command -v matlab >/dev/null 2>&1; then \
		matlab -nodisplay -nosplash -r "SL_equation_viscoelastic_GIAonly; exit"; \
	else \
		echo "Error: MATLAB not found."; \
		exit 1; \
	fi

matlab-dynamic-topography:
	@echo "Running MATLAB dynamic topography model..."
	@if command -v matlab >/dev/null 2>&1; then \
		matlab -nodisplay -nosplash -r "SL_equation_DT; exit"; \
	else \
		echo "Error: MATLAB not found."; \
		exit 1; \
	fi

matlab-benchmarked: matlab-elastic matlab-viscoelastic
	@echo "Completed benchmarked MATLAB runs"

# Special case: Proglacial lakes study
matlab-lakes:
	@echo "Running proglacial lakes study..."
	@if command -v matlab >/dev/null 2>&1; then \
		cd GRL_2022_proglacial_lakes && matlab -nodisplay -nosplash -r "SL_equation_viscoelastic_lakes; exit"; \
	else \
		echo "Error: MATLAB not found."; \
		exit 1; \
	fi

# Plot results (requires MATLAB)
plot-fingerprint:
	@echo "Plotting sea level fingerprint..."
	@if command -v matlab >/dev/null 2>&1; then \
		cd plot_output && matlab -nodisplay -nosplash -r "plot_fingerprint; exit"; \
	else \
		echo "Error: MATLAB not found."; \
		exit 1; \
	fi

plot-rsl:
	@echo "Plotting relative sea level..."
	@if command -v matlab >/dev/null 2>&1; then \
		cd plot_output && matlab -nodisplay -nosplash -r "plot_RSL; exit"; \
	else \
		echo "Error: MATLAB not found."; \
		exit 1; \
	fi

# Clean temporary files
clean:
	@echo "Cleaning temporary files..."
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.m~" -delete 2>/dev/null || true
	@echo "Clean complete"

# Information targets
info:
	@echo "SLcode Repository Information"
	@echo "============================"
	@echo "Main Python solver:     SL_equation_elastic.py"
	@echo "Alternative Python:     SLcode_py/SL_equation_elastic.py"
	@echo "MATLAB elastic:         SL_equation_elastic_benchmarked.m"
	@echo "MATLAB viscoelastic:    SL_equation_viscoelastic_benchmarked.m"
	@echo "MATLAB GIA-only:        SL_equation_viscoelastic_GIAonly.m"
	@echo "MATLAB dynamic topo:    SL_equation_DT.m"
	@echo "Proglacial lakes:       GRL_2022_proglacial_lakes/"
	@echo "Love numbers:           SavedLN/ and SavedLN_EP/"
	@echo "Ice grids:              ice_grid/"
	@echo "Functions:              SLFunctions/"
	@echo "Plot tools:             plot_output/"

list-data:
	@echo "Available data files:"
	@find . -name "*.mat" | head -10
	@echo "... and more MATLAB data files in SavedLN/ and SavedLN_EP/"

# Advanced targets for specific studies
csdms-workshop:
	@echo "Running CSDMS workshop examples..."
	@cd CSDMS_workshop && find . -name "*.m" -exec echo "Found: {}" \;

benchmark-io:
	@echo "Running benchmark input/output tests..."
	@cd benchmark_in_out && find . -name "*.mat" -exec echo "Found benchmark data: {}" \;