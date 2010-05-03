CC=gcc
#CPPFLAGS=
CFLAGS=-fPIC -Wall -O2
LDFLAGS=-shared
LDLIBS=-lgsl -lgslcblas

all: libar.so

libar.so: libar.o:
	$(CC) $(LDFLAGS) $(LDLIBS) -Wl,-soname,$@ -o $@ $<

libar.o: libar.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.so


# install/%.la: %.la
#         libtool --mode=install \
#         install -c $(notdir $@) $(libdir)/$(notdir $@)
# install: $(addprefix install/,$(LIBS))
#         libtool --mode=finish $(libdir)
