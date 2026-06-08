CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -Wno-unused-parameter
TARGET   = bin/solver
TARGER_DIR = bin
SRCDIR   = src

SRC  = $(SRCDIR)/main.cpp
HDRS = $(SRCDIR)/config.h $(SRCDIR)/types.h $(SRCDIR)/parser.h $(SRCDIR)/bstree.h \
       $(SRCDIR)/ft_estimator.h $(SRCDIR)/channel.h $(SRCDIR)/router.h \
       $(SRCDIR)/floorplan.h $(SRCDIR)/sa_optimizer.h $(SRCDIR)/output.h

$(TARGET): $(SRC) $(HDRS)
	@mkdir -p $(TARGER_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: clean
