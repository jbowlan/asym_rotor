CC=gcc
CPPFLAGS=
CFLAGS=-Wall -O3
LDFLAGS=
LDLIBS=gsl gslcblas

define compile_rule
        libtool --mode=compile \
        $(CC) $(CFLAGS) $(CPPFLAGS) -c $<
endef
define link_rule
        libtool --mode=link \
        $(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)
endef

LIBS = libar.la
libar_OBJS = asym_rotor.lo

%.lo: %.c
        $(call compile_rule)

libar.la: $(libmystuff_OBJS)
        $(call link_rule)

# install/%.la: %.la
#         libtool --mode=install \
#         install -c $(notdir $@) $(libdir)/$(notdir $@)
# install: $(addprefix install/,$(LIBS))
#         libtool --mode=finish $(libdir)
