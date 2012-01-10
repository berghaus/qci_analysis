#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT(analysis, 1.0, frank at uvic.ca)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lCint':
AC_CHECK_LIB([Cint], [main])
# FIXME: Replace `main' with a function in `-lCore':
AC_CHECK_LIB([Core], [main])
# FIXME: Replace `main' with a function in `-lGpad':
AC_CHECK_LIB([Gpad], [main])
# FIXME: Replace `main' with a function in `-lGraf':
AC_CHECK_LIB([Graf], [main])
# FIXME: Replace `main' with a function in `-lGraf3d':
AC_CHECK_LIB([Graf3d], [main])
# FIXME: Replace `main' with a function in `-lGui':
AC_CHECK_LIB([Gui], [main])
# FIXME: Replace `main' with a function in `-lHist':
AC_CHECK_LIB([Hist], [main])
# FIXME: Replace `main' with a function in `-lMathCore':
AC_CHECK_LIB([MathCore], [main])
# FIXME: Replace `main' with a function in `-lMatrix':
AC_CHECK_LIB([Matrix], [main])
# FIXME: Replace `main' with a function in `-lMinuit2':
AC_CHECK_LIB([Minuit2], [main])
# FIXME: Replace `main' with a function in `-lNet':
AC_CHECK_LIB([Net], [main])
# FIXME: Replace `main' with a function in `-lPhysics':
AC_CHECK_LIB([Physics], [main])
# FIXME: Replace `main' with a function in `-lPostscript':
AC_CHECK_LIB([Postscript], [main])
# FIXME: Replace `main' with a function in `-lRIO':
AC_CHECK_LIB([RIO], [main])
# FIXME: Replace `main' with a function in `-lRint':
AC_CHECK_LIB([Rint], [main])
# FIXME: Replace `main' with a function in `-lThread':
AC_CHECK_LIB([Thread], [main])
# FIXME: Replace `main' with a function in `-lTree':
AC_CHECK_LIB([Tree], [main])
# FIXME: Replace `main' with a function in `-ldl':
AC_CHECK_LIB([dl], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
