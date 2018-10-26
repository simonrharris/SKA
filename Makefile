CXX = g++
CXXFLAGS = -g -std=c++0x -O2

RM = rm -f
MKDIR = mkdir
ECHO = echo
CP = cp
INSTALL = install
dir_guard=@mkdir -p $(@D)

PREFIX = /usr/local
SOURCEDIR = src
HEADERDIR = src
BUILDDIR = build
BINARYDIR = bin

SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))

LDFLAGS = -lz

# the build target executable:
BINARY = ska

all: $(BINARYDIR)/$(BINARY)

$(BINARYDIR)/$(BINARY): $(OBJECTS)
	$(dir_guard)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SOURCEDIR)/%.cpp
	$(dir_guard)
	$(CXX) $(CXXFLAGS) -I$(HEADERDIR) -I$(SOURCEDIR) -c $< -o $@

.phony: install
install: all
	$(INSTALL) -d $(PREFIX)/bin
	$(INSTALL) $(BINARYDIR)/$(BINARY) $(PREFIX)/bin

.phony: clean
clean:
	$(RM) $(OBJECTS)

.phony: distclean
distclean: clean
	$(RM) $(BINARYDIR)/$(BINARY)

help:
	@$(ECHO) "Targets:"
	@$(ECHO) "all     - build and compile what is necessary"
	@$(ECHO) "install   - install binary"
	@$(ECHO) "clean   - cleanup old .o files"
	@$(ECHO) "distclean   - cleanup old binary"
