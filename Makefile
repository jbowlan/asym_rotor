CC=gcc
#CPPFLAGS=
CFLAGS=-fPIC -Wall -O2 -I/opt/local/include
LDFLAGS=-shared -L/opt/local/lib
LDLIBS=-lgsl -lgslcblas

all: libar.so

libar.so: libar.o
	$(CC) $(LDFLAGS) $(LDLIBS) -Wl -o libar.dylib libar.o

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.dylib


# install/%.la: %.la
#         libtool --mode=install \
#         install -c $(notdir $@) $(libdir)/$(notdir $@)
# install: $(addprefix install/,$(LIBS))
#         libtool --mode=finish $(libdir)
