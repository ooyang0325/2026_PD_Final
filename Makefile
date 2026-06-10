CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -Wno-unused-parameter
TARGET   = bin/solver
TARGER_DIR = bin
SRCDIR   = src
VENDDIR  = $(SRCDIR)/vendor
OBJDIR   = build

SRC  = $(SRCDIR)/main.cpp
HDRS = $(SRCDIR)/config.h $(SRCDIR)/types.h $(SRCDIR)/parser.h $(SRCDIR)/bstree.h \
       $(SRCDIR)/ft_estimator.h $(SRCDIR)/channel.h $(SRCDIR)/router.h \
       $(SRCDIR)/floorplan.h $(SRCDIR)/sa_optimizer.h $(SRCDIR)/legalize_loop.h \
       $(SRCDIR)/analytical_legalizer.h $(SRCDIR)/output.h \
       $(SRCDIR)/mp_optimizer.h $(VENDDIR)/fft.h

OBJS = $(OBJDIR)/main.o \
       $(OBJDIR)/fft.o \
       $(OBJDIR)/fftsg.o \
       $(OBJDIR)/fftsg2d.o

$(TARGET): $(OBJS)
	@mkdir -p $(TARGER_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(OBJDIR)/main.o: $(SRC) $(HDRS)
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $@

$(OBJDIR)/fft.o: $(VENDDIR)/fft.cpp $(VENDDIR)/fft.h
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $(VENDDIR)/fft.cpp -o $@

$(OBJDIR)/fftsg.o: $(VENDDIR)/fftsg.cpp $(VENDDIR)/fft.h
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $(VENDDIR)/fftsg.cpp -o $@

$(OBJDIR)/fftsg2d.o: $(VENDDIR)/fftsg2d.cpp $(VENDDIR)/fft.h
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $(VENDDIR)/fftsg2d.cpp -o $@

clean:
	rm -f $(TARGET)
	rm -rf $(OBJDIR)

.PHONY: clean
