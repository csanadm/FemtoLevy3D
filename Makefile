ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit2
ROOTGLIBS    := $(shell root-config --glibs) -lMinuit2

WrkDir := $(shell pwd)
SrcDir        = ./
ObjDir        = ./
DepDir        = ./deps
ExeDir        = ./

PROGRAMS = $(notdir $(basename $(shell ls $(SrcDir)/*.cc)))

CXX           = g++
CFLAGS        = -O3 -Wall -fPIC -fno-inline -std=c++11 -I$(WrkDir)/$(SrcDir) -I/star/u/kincses/boost/boost_1_82_0/
#CFLAGS        = -O3 -Wall -fPIC -fno-inline -std=c++14 -I$(WrkDir)/$(SrcDir) -I/star/u/kincses/boost/boost_1_82_0/
LD            = g++
LDFLAGS       = -O3 
RFLAGS        = $(CFLAGS) $(ROOTCFLAGS) -I/star/u/kincses/boost/boost_1_82_0/
RLIBS         = $(ROOTLIBS) $(SYSLIBS)
RGLIBS        = $(ROOTGLIBS) $(SYSLIBS)

COMMON_SOURCES = $(shell ls $(SrcDir)/*.cpp)
COMMON_OBJECTS = $(addprefix $(ObjDir)/, $(addsuffix .o,$(notdir $(basename $(COMMON_SOURCES)))))
SOURCES = $(addprefix $(SrcDir)/,$(addsuffix .cc,$(PROGRAMS))) $(COMMON_SOURCES)
ALL_SOURCES = $(sort $(SOURCES))

define COMPILE_TEMPLATE
-include $(DepDir)/$(notdir $(basename $(1))).d
$(ObjDir)/$(notdir $(basename $(1))).o:
	@echo ""
	$(CXX) $(RFLAGS) -c -MD -MP -MF $(DepDir)/$(notdir $(basename $(1))).d $(1) -o $$@
endef
$(foreach source, $(ALL_SOURCES), $(eval $(call COMPILE_TEMPLATE,$(source))))

define FINAL_LINK_TEMPLATE
$(ExeDir)/$(1).exe: $(2) $(ObjDir)/$(1).o
	@echo ""
	$(LD) $(LDFLAGS) $$^ $(RGLIBS) -o $$@
endef
$(foreach program, $(PROGRAMS), $(eval $(call FINAL_LINK_TEMPLATE,$(program),$(COMMON_OBJECTS))))

all: $(addprefix $(ExeDir)/,$(addsuffix .exe, $(PROGRAMS)))

clean:
	@rm -f $(ExeDir)/*.exe $(ObjDir)/*.o $(DepDir)/*.d
