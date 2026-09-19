PGR_SRC       += $(wildcard *.cc)

PGR_OBJS   = $(patsubst %.cc, %.o, $(PGR_SRC))

PROGRAMS   = $(patsubst %.o,%, $(PGR_OBJS))

CXXFLAGS += `root-config --cflags` -I.
LDFLAGS  += `root-config --ldflags`
LDLIBS  += `root-config --glibs`

all: $(PROGRAMS) $(OBJS) 

clean:
	@rm -f $(OBJS) $(PGR_OBJS) $(PROGRAMS) *~ *Dict.* $(LIB)

$(PROGRAMS): %: %.o $(OBJS) 
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@
%.o:%.cxx
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o:%.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o:%.c
	$(CC) -c $(CCFLAGS) $< -o $@

%Dict.cc:%.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

