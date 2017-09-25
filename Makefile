#------------------------------------------------------------------------------#
#                     Makefile for ED and Miami-LU models                      #
#------------------------------------------------------------------------------#


CXX = g++

INC = -I/opt/local/include -I/opt/local/include/db53
#INC = -I$(DBXML_INCLUDE) -I/lustre/data/fisk/local/include -I$(NETCDFINC)

LIB = -L/opt/local/lib -L/opt/local/lib/db53
#LIB = -L$(DBXML_LIB) -L/lustre/data/fisk/local/lib -L$(NETCDFLIB)

# add  -g for debug, -g3 for debug and optomize, -O2 for optimize
#      -Wall for full warnings
#      -pg for profiling
CXXFLAGS = $(INC) -c -Wall -g
LDFLAGS = $(LIB) -lm -lnetcdf -lnetcdf_c++ -ltbb -ldb_cxx -lconfig++

CMN_SRCS = site.cc patch.cc miami.cc belowgrnd.cc \
           disturbance.cc fire.cc landuse.cc read_site_data.cc init_data.cc \
           outputter.cc print_output.cc restart.cc readconfiguration.cc

EDM_SRCS = cohort.cc growth.cc allometry.cc phenology.cc mortality.cc \
           mechanism.cc odeint.cc


# Can't use target specific variables because they aren't available 
# when needed for dependency includes
ifeq ($(MAKECMDGOALS),mlu)
	TGT = mlu
   CXXFLAGS += -DMIAMI_LU -DMAIN
	SRCS = $(CMN_SRCS) main.cc
else ifeq ($(MAKECMDGOALS),ed_mpi)
   TGT = ed_mpi
   CXX = mpicxx
   CXXFLAGS += -DED -DUSEMPI -DMAIN
   SRCS = $(CMN_SRCS) $(EDM_SRCS) edmpi.cc main.cc
else ifeq ($(MAKECMDGOALS),libmlu.a)
	TGT = libmlu
   CXXFLAGS += -DMIAMI_LU -DCOUPLED
   SRCS = $(CMN_SRCS) edm_ied_interface.cc
else ifeq ($(MAKECMDGOALS),libed.a)
	TGT = libed
   CXXFLAGS += -DED -DCOUPLED
   SRCS = $(CMN_SRCS) $(EDM_SRCS) edm_ied_interface.cc
else
	TGT = edlu
   CXXFLAGS += -DED -DMAIN
   SRCS = $(CMN_SRCS) $(EDM_SRCS) main.cc
endif

OBJS = $(SRCS:.cc=.o)

all: edlu

edlu ed_mpi mlu: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@

libed.a libmlu.a: $(OBJS)
	ar rcs $@ $(OBJS)

%.tgt: 
	@rm -f *.tgt; \
	touch $@; 

# make sure that .o and .d files get rebuilt if we switch targets
$(OBJS) $(SRCS:.cc=.d): $(TGT).tgt

# automatically determine header dependencies 
%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
   -include $(SRCS:.cc=.d)
endif


.PHONY: clean
clean:
	- rm -f *.d *.o core.* edlu mlu

