DIR_SRC = ./src
DIR_OBJ = ./obj
BINDIR = /usr/local/bin
HTSLIB_INC ?=
HTSLIB_LIB ?=

SRC = $(wildcard $(DIR_SRC)/*.cpp)
OBJ = $(patsubst %.cpp,$(DIR_OBJ)/%.o,$(notdir $(SRC)))

TARGET = gencore
BIN_TARGET = $(TARGET)

CXX ?= g++
CPPFLAGS += -I$(DIR_SRC)
ifneq ($(wildcard ./inc),)
CPPFLAGS += -I./inc
endif
ifneq ($(HTSLIB_INC),)
CPPFLAGS += -I$(HTSLIB_INC)
endif
ifneq ($(HTSLIB_LIB),)
LDFLAGS += -L$(HTSLIB_LIB)
endif
CXXFLAGS += -std=c++11 -g -O3
LDLIBS += -lhts -lz -pthread

$(BIN_TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $@

$(DIR_OBJ)/%.o: $(DIR_SRC)/%.cpp make_obj_dir
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(DIR_OBJ)/*.o
	rm -f $(TARGET)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."
