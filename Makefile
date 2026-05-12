CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -Wno-unused-parameter
TARGET   = solver
SRCDIR   = src

SRC  = $(SRCDIR)/main.cpp
HDRS = $(SRCDIR)/types.h $(SRCDIR)/parser.h $(SRCDIR)/sequence_pair.h \
       $(SRCDIR)/channel.h $(SRCDIR)/router.h $(SRCDIR)/floorplan.h \
       $(SRCDIR)/sa_optimizer.h $(SRCDIR)/output.h

$(TARGET): $(SRC) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: clean
