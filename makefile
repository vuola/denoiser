# denoiser - build, test and remote-test targets.
#
#   make               build the program            -> build/denoiser
#   make test          build and run the unit tests -> build/run_tests
#   make reference     build the original, unchanged filter -> build/denoiser_reference
#                      and its one-worker variant          -> build/denoiser_reference_1w
#   make clean         remove all build output
#   make remote-test   run the tests of the last commit on the remote machine
#   make help          list the targets
#
# Options (append to the command, e.g. "make test SANITIZE=1"):
#   SANITIZE=1   build with address and undefined-behaviour checks into
#                build-san/ (Linux only; MinGW on Windows has no sanitizers)
#   REMOTE=name  SSH host alias of the remote machine (default: gpu)
#
# On Windows make needs the Unix tools of Git for Windows (sh, mkdir, rm) on the
# PATH; it works from Git Bash and from PowerShell.

CXX      := g++
# No -march=native and no -ffast-math: results must be comparable between
# machines (README section 1.4).
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -Iinclude
OPENMP   := -fopenmp
BUILD    := build

ifeq ($(SANITIZE),1)
  CXXFLAGS += -g -fno-omit-frame-pointer -fsanitize=address,undefined
  BUILD    := build-san
endif

# Windows: programs end in .exe, and the compiler's runtime libraries are
# linked into the program (-static). Otherwise the program picks up whichever
# libstdc++ DLL comes first on the PATH; Git for Windows ships an older one,
# and the program then fails to start.
ifeq ($(OS),Windows_NT)
  EXE     := .exe
  LDFLAGS := -static
else
  EXE     :=
  LDFLAGS :=
endif

# Remote machine for "make remote-test": SSH host alias and the clone's folder
# there. The alias is defined in ~/.ssh/config, outside the repository.
REMOTE     := gpu
REMOTE_DIR := denoiser

# Library code: everything in src/ except the command-line entry point.
LIB_SRC  := $(filter-out src/main.cpp,$(wildcard src/*.cpp))
TEST_SRC := $(wildcard tests/*.cpp)

LIB_OBJ  := $(LIB_SRC:%.cpp=$(BUILD)/%.o)
MAIN_OBJ := $(BUILD)/src/main.o
TEST_OBJ := $(TEST_SRC:%.cpp=$(BUILD)/%.o)

PROGRAM  := $(BUILD)/denoiser$(EXE)
TESTS    := $(BUILD)/run_tests$(EXE)

.PHONY: all test reference clean remote-test help

all: $(PROGRAM)

$(PROGRAM): $(MAIN_OBJ) $(LIB_OBJ)
	$(CXX) $(CXXFLAGS) $(OPENMP) $^ $(LDFLAGS) -o $@

$(TESTS): $(TEST_OBJ) $(LIB_OBJ)
	$(CXX) $(CXXFLAGS) $(OPENMP) $^ $(LDFLAGS) -o $@

# Tests run from the repository root so that they find tests/data/.
test: $(TESTS)
	./$(TESTS)

# One rule for every object file. -MMD -MP writes a .d file listing the headers
# each object depends on, so a changed header rebuilds what uses it.
$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(OPENMP) -Ithird_party -MMD -MP -c $< -o $@

-include $(LIB_OBJ:.o=.d) $(MAIN_OBJ:.o=.d) $(TEST_OBJ:.o=.d)

# The original filter, compiled without any change to its source. It needs
# C++20 and uses std::thread. The source is local reference material that is
# not in git (README section 1.3); the target fails if the file is absent.
#
# Two programs are built from the same source:
#   denoiser_reference     as the author built it: (processors - 1) workers;
#                          the result differs slightly from run to run
#   denoiser_reference_1w  one worker (tests/reference/one_worker.cpp replaces
#                          the processor count); the result is repeatable
REFERENCE_SRC := ORIGINAL_CODE/denoiserVersionA.cpp
REFERENCE     := $(BUILD)/denoiser_reference$(EXE)
REFERENCE_1W  := $(BUILD)/denoiser_reference_1w$(EXE)
REFERENCE_CXX := $(CXX) -O2 -std=c++20

# With the static runtime on Windows the linker must be told that the
# replacement of the processor count is intended.
ifeq ($(OS),Windows_NT)
  ONE_WORKER_LDFLAGS := -Wl,--allow-multiple-definition
else
  ONE_WORKER_LDFLAGS :=
endif

reference: $(REFERENCE) $(REFERENCE_1W)

$(REFERENCE): $(REFERENCE_SRC)
	@mkdir -p $(BUILD)
	$(REFERENCE_CXX) $^ $(LDFLAGS) -pthread -o $@

$(REFERENCE_1W): $(REFERENCE_SRC) tests/reference/one_worker.cpp
	@mkdir -p $(BUILD)
	$(REFERENCE_CXX) $^ $(LDFLAGS) $(ONE_WORKER_LDFLAGS) -pthread -o $@

clean:
	rm -rf build build-san

# Sends the last commit (HEAD) to the remote clone as branch "remote-test",
# checks it out there and runs the tests. Uncommitted changes are not sent.
# Nothing is pushed to GitHub.
remote-test:
	git push --force $(REMOTE) HEAD:refs/heads/remote-test
	ssh $(REMOTE) "cd $(REMOTE_DIR) && git checkout --quiet --force --detach remote-test && make test"

help:
	@sed -n '2,14p' makefile
