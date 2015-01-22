CC = g++
CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++11
# this is very unportable
LDFLAGS := $(LDFLAGS) $(shell pkg-config --libs-only-L jellyfish-2.0) -L$(BOOST_ROOT)/lib 
LDLIBS = -lboost_timer -lboost_chrono -lboost_system -lboost_program_options  $(shell pkg-config --libs-only-l jellyfish-2.0)

all: bin/jellyfish-mirror
bin/jellyfish-mirror: jellyfish-mirror
	cp $< $@

jellyfish-mirror: jellyfish-mirror.cc
clean:
	rm -f bin/* jellyfish-mirror
