#https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

# CC  = clang++-12
# CXX = clang++-12
# LD  = clang++-12

CC  = g++
CXX = g++
LD  = g++

RM = rm -rf

OPTIM_FLAGS = -O3 -lm
GTEST_FLAGS = -lgtest -lgtest_main -pthread
DEBUG_FLAGS = -g
# SANITISE_FLAGS = -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -static-libasan
OTHER_FLAGS = -std=c++17 -MD

COMMON_FLAGS = $(OPTIM_FLAGS) $(OTHER_FLAGS) $(DEBUG_FLAGS) $(SANITISE_FLAGS)
CFLAGS   = $(COMMON_FLAGS)
CXXFLAGS = $(COMMON_FLAGS)
LDFLAGS  = $(COMMON_FLAGS)

OBJS = main.o BinArray.o
TEST_DIR = tests
SANDBOX_DIR = sandbox

all: main.out test.out gen_test.out

# Main -------------------------
main.out: main.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

# Test -------------------------
test.out: test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(GTEST_FLAGS)

gen_test.out: gen_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

full_solve.out: full_solve.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean -------------------------
tmp_clean:
	$(RM) *.o *.d

clean: tmp_clean
	$(RM) *.out

clean_test:
	$(RM) test_ref.txt test.txt

-include main.d test.d gen_test.d full_solve.d

.PHONY: all clean_temp clean clean_test