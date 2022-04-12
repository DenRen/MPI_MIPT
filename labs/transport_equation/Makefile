#https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html

CC  = mpic++
CXX = mpic++
LD  = mpic++

RM = rm -rf

OPTIM_FLAGS = -O3 -lm -DNDEBUG
GTEST_FLAGS = -lgtest -lgtest_main -pthread
# DEBUG_FLAGS = -g3
# SANITISE_FLAGS = -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -static-libasan
OTHER_FLAGS = -std=c++17 -MD -pthread

COMMON_FLAGS = $(OPTIM_FLAGS) $(OTHER_FLAGS) $(DEBUG_FLAGS) $(SANITISE_FLAGS)
CFLAGS   = $(COMMON_FLAGS)
CXXFLAGS = $(COMMON_FLAGS)
LDFLAGS  = $(COMMON_FLAGS)

OBJS = treq_lib.o

all: main.out test.out

# Main -------------------------
main.out: main.o $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Test -------------------------
test.out: test.o $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(GTEST_FLAGS)

# Clean -------------------------
tmp_clean:
	$(RM) *.o *.d *.gch

clean: tmp_clean
	$(RM) *.out

clean_test:
	$(RM) test_ref.txt test.txt

-include treq_lib.d print_lib.d func_lib.d copy_transform.d

.PHONY: all clean_temp clean clean_test