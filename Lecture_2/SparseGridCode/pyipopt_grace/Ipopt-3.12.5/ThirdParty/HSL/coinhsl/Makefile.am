AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libcoinhsl.la
include_HEADERS = CoinHslConfig.h

pkgconfigdir 	= $(libdir)/pkgconfig
pkgconfig_DATA = coinhsl.pc

libcoinhsl_la_SOURCES = common/deps90.f90 \
								common/deps.f \
								common/dump.f90

EXTRA_DIST = LICENCE

if HAVE_MC19
libcoinhsl_la_SOURCES += mc19/mc19d.f
endif

if HAVE_MA27
libcoinhsl_la_SOURCES += ma27/ma27d.f
endif

if HAVE_MA28
libcoinhsl_la_SOURCES += ma28/ma28d.f
endif

if HAVE_MA57
libcoinhsl_la_SOURCES += ma57/ma57d.f
endif

if HAVE_MA77
include_HEADERS 		 += hsl_ma77/C/hsl_ma77d.h
libcoinhsl_la_SOURCES += hsl_ma77/hsl_ma77d.f90 \
								 hsl_ma77/C/hsl_ma77d_ciface.f90
endif

if HAVE_MA86
include_HEADERS       += hsl_ma86/C/hsl_ma86d.h
libcoinhsl_la_SOURCES += hsl_ma86/hsl_ma86d.f90 \
								 hsl_ma86/C/hsl_ma86d_ciface.f90
endif

if HAVE_MC68
include_HEADERS 		 += hsl_mc68/C/hsl_mc68i.h
#libcoinhsl_la_SOURCES += hsl_mc68/hsl_mc68i.f90 \
#								 hsl_mc68/C/hsl_mc68i_ciface.f90
libcoinhsl_la_SOURCES += hsl_mc68/C/hsl_mc68i_ciface.f90
endif

if HAVE_MA97
include_HEADERS 		 += hsl_ma97/C/hsl_ma97d.h
libcoinhsl_la_SOURCES += hsl_ma97/hsl_ma97d.f90 \
								 hsl_ma97/C/hsl_ma97d_ciface.f90
endif

# Automake doesn't know much about modules, so get rid of them ourself
clean-local:
	-rm *.mod

hsl_ma86d_ciface.$(OBJEXT): hsl_ma86d.$(OBJEXT)
hsl_ma97d_ciface.$(OBJEXT): hsl_ma97d.$(OBJEXT)
hsl_mc68i_ciface.$(OBJEXT): hsl_mc68i.$(OBJEXT)
