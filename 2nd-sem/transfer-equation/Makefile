
.PHONY: all time dump  

GENERAL_OPTS=-std=c++2a -g -I matrix/ 
TIME_OPTS=-D"TIME"
DUMP_OPTS=-D"DUMP"

COMPILE_SEQUENTIAL    = g++ sequential.cpp helper.cpp -o sequential.out
COMPILE_PARALLEL = mpic++ parallel.cpp helper.cpp -o parallel.out 

all: 
	$(COMPILE_SEQUENTIAL) $(GENERAL_OPTS) 
	$(COMPILE_PARALLEL) $(GENERAL_OPTS) 

time: 
	$(COMPILE_SEQUENTIAL) $(GENERAL_OPTS) $(TIME_OPTS)
	$(COMPILE_PARALLEL) $(GENERAL_OPTS) $(TIME_OPTS)

dump:
	$(COMPILE_SEQUENTIAL) $(GENERAL_OPTS) $(DUMP_OPTS) 
	$(COMPILE_PARALLEL) $(GENERAL_OPTS) $(DUMP_OPTS)
