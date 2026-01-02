# Makefile wrapper for CMake build system
# Provides convenient "make" interface while using CMake underneath
export CMAKE_EXPORT_COMPILE_COMMANDS := ON

.PHONY: all clean debug help

# Check if CMake is available
CMAKE := $(shell command -v cmake 2> /dev/null)

# Allow forcing root-config mode via: make ROOTCONFIG=1
# or by setting FORCE_ROOT_CONFIG environment variable
ifdef ROOTCONFIG
  CMAKE_ROOT_FLAGS := -DFORCE_ROOT_CONFIG=ON
else ifdef FORCE_ROOT_CONFIG
  CMAKE_ROOT_FLAGS := -DFORCE_ROOT_CONFIG=ON
else
  CMAKE_ROOT_FLAGS :=
endif

# Extract parallel job count from MAKEFLAGS if present
# This allows "make -j8" to pass through to cmake --build
JOBS := $(patsubst -j%,%,$(filter -j%,$(MAKEFLAGS)))
ifneq ($(JOBS),)
  CMAKE_BUILD_FLAGS := --parallel $(JOBS)
else
  # No -j flag specified, let CMake decide (usually # of cores)
  CMAKE_BUILD_FLAGS := --parallel
endif

# Default target: build the project
all:
ifndef CMAKE
	$(error CMake is not installed or not in PATH. Please install CMake 3.16+ and try again. See https://cmake.org/download/)
endif
	@test -d build || cmake -B build $(CMAKE_ROOT_FLAGS)
	@cmake --build build $(CMAKE_BUILD_FLAGS)

# Build with debug symbols and no optimization
debug:
ifndef CMAKE
	$(error CMake is not installed or not in PATH. Please install CMake 3.16+ and try again. See https://cmake.org/download/)
endif
	@cmake -B build -DCMAKE_BUILD_TYPE=Debug $(CMAKE_ROOT_FLAGS)
	@cmake --build build $(CMAKE_BUILD_FLAGS)

# Force reconfiguration and rebuild
reconfigure:
ifndef CMAKE
	$(error CMake is not installed or not in PATH. Please install CMake 3.16+ and try again. See https://cmake.org/download/)
endif
	@cmake -B build $(CMAKE_ROOT_FLAGS)
	@cmake --build build $(CMAKE_BUILD_FLAGS)

# Clean build artifacts
clean:
	@rm -rf build

# Show help message
help:
	@echo "XSecAnalyzer Build Targets:"
	@echo "  make                - Build project (default, auto-configures if needed)"
	@echo "  make -j8            - Build with 8 parallel jobs (or any number)"
	@echo "  make ROOTCONFIG=1   - Force use of root-config instead of find_package"
	@echo "  make debug          - Build with debug symbols"
	@echo "  make reconfigure    - Force CMake reconfiguration"
	@echo "  make clean          - Remove build directory"
	@echo "  make help           - Show this message"
	@echo ""
	@echo "You can also set FORCE_ROOT_CONFIG=1 as an environment variable"
