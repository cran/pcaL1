AC_INIT([pcaL1], [1.5.4], [jpbrooks@vcu.edu])



dnl # --------------------------------------------------------------------
dnl # global (environment) variables
dnl # --------------------------------------------------------------------
dnl
dnl # PKG_CPPFLAGS       C preprocessor flags -I...
dnl # PKG_CFLAGS         C compiler flags
dnl # PKG_LIBS           libraries -l...


dnl # --------------------------------------------------------------------
dnl # R compiler and flags
dnl # --------------------------------------------------------------------

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "could not determine R_HOME"
    exit 1
fi

dnl # GNU R default settings
CC=`"${R_HOME}/bin/R" CMD config CC`
#CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP



dnl # --------------------------------------------------------------------
dnl # test for global variables
dnl # --------------------------------------------------------------------

ARGS_OK="FALSE"

if test [ -n "$PKG_CPPFLAGS" -a -n "$PKG_LIBS" ] ; then
    ARGS_OK="TRUE"
else
    if test [ -n "$PKG_CPPFLAGS" -o -n "$PKG_LIBS" ] ; then
        AC_MSG_ERROR([Please set both variables: PKG_CPPFLAGS AND PKG_LIBS])
    fi
fi


dnl # --------------------------------------------------------------------
dnl # get arguments
dnl # --------------------------------------------------------------------

AC_ARG_WITH(clp-include,
            AC_HELP_STRING([--with-clp-include=PATH],
                           [location of COIN-OR Clp include directory [[/usr/local/include]]]
                          ),
            [CLP_INCLUDE=$withval],
            [CLP_INCLUDE=""])

AC_ARG_WITH(clp-lib,
            AC_HELP_STRING([--with-clp-lib=PATH],
                           [location of COIN-OR Clp callable library [[usr/local/lib]]]
                          ),
            [CLP_LIB=$withval],
            [CLP_LIB=""])


dnl # --------------------------------------------------------------------
dnl # pkg-config
dnl # --------------------------------------------------------------------

dnl # include and library directories
if test [ "$ARGS_OK" = "FALSE" ] ; then
    AC_PATH_PROG([PKG_CONFIG], [pkg-config], [], [$PATH:/usr/local/bin])

    if test [ "${PKG_CONFIG}" != "" ] ; then
        if "${PKG_CONFIG}" --exists clp; then
    
            CLP_LIBS=`${PKG_CONFIG} --libs clp`    
            CLP_CFLAGS=`${PKG_CONFIG} --cflags clp`    
    
            if test [ -z "$CLP_LIBS" -a -z "$CLP_CFLAGS" ] ; then
                AC_MSG_NOTICE([can not get information about clp from pkg-config])
            else
                AC_MSG_NOTICE([set compile and link flags according to pkg-config])
                AC_MSG_NOTICE([compile flags are $CLP_CFLAGS])
                AC_MSG_NOTICE([link flags are $CLP_LIBS])
                PKG_LIBS="${CLP_LIBS}"
                PKG_CPPFLAGS="${CLP_CFLAGS}"
                ARGS_OK="TRUE"
            fi
        else
            AC_MSG_NOTICE([pkg-config file for clp is not available])
        fi
    fi
fi


dnl # --------------------------------------------------------------------
dnl # test arguments
dnl # --------------------------------------------------------------------

dnl # include and library directories

dnl # include and library directories
if test [ "$ARGS_OK" = "FALSE" ] ; then

    if test [ -n "$CLP_INCLUDE" -a -n "$CLP_LIB" ] ; then
    
        if test [ -d "$CLP_INCLUDE" ] ; then
            if test [ -d "$CLP_INCLUDE/coin" ] ; then
                PKG_CPPFLAGS="-I${CLP_INCLUDE} -I${CLP_INCLUDE}/coin"
            else
                PKG_CPPFLAGS="-I${CLP_INCLUDE}"
            fi
        else
            AC_MSG_ERROR([directory $CLP_INCLUDE does not exist])
        fi
    
        if test [ -d "$CLP_LIB" ] ; then
            if test [ -d "$CLP_LIB/coin" ] ; then
                PKG_LIBS="-L${CLP_LIB} -L${CLP_LIB}/coin"
            else
                PKG_LIBS="-L${CLP_LIB}"
            fi
        else
            AC_MSG_ERROR([directory $CLP_LIB does not exist])
        fi
    
    else
    
    dnl # if no arguments are given, try to find the clp executable. If is not
    dnl # found, use /usr/local
    
        AC_PREFIX_PROGRAM(clp)
    
        if test [ "$prefix" = "NONE" ] ; then
    
            AC_MSG_NOTICE([No prefix given, trying /usr and /usr/local])
    
            if test [ -d "/usr/include" -o -d "/usr/local/include" ] ; then
                PKG_CPPFLAGS="-I/usr/include -I/usr/local/include -I/usr/include/coin -I/usr/local/include/coin"
            else
                AC_MSG_ERROR([directories /usr/include and /usr/local/include do not exist])
            fi

            if test [ -d "/usr/local/lib" -o -d "/usr/local/lib64" -o -d "/usr/lib" -o -d "/usr/lib64"] ; then
                PKG_LIBS=""
                if test [ -d "/usr/local/lib" ] ; then
                    PKG_LIBS="$PKG_LIBS -L/usr/local/lib"
                    if test [ -d "/usr/local/lib/coin" ] ; then
                        PKG_LIBS="$PKG_LIBS -L/usr/local/lib/coin"
                    fi
                fi
                if test [ -d "/usr/local/lib64" ] ; then
                    PKG_LIBS="$PKG_LIBS -L/usr/local/lib64"
                    if test [ -d "/usr/local/lib64/coin" ] ; then
                        PKG_LIBS="$PKG_LIBS -L/usr/local/lib64/coin"
                    fi
                fi
                if test [ -d "/usr/lib" ] ; then
                    PKG_LIBS="$PKG_LIBS -L/usr/lib"
                    if test [ -d "/usr/lib/coin" ] ; then
                        PKG_LIBS="$PKG_LIBS -L/usr/lib/coin"
                    fi
                fi
                if test [ -d "/usr/lib64" ] ; then
                    PKG_LIBS="$PKG_LIBS -L/usr/lib64"
                    if test [ -d "/usr/lib64/coin" ] ; then
                        PKG_LIBS="$PKG_LIBS -L/usr/lib64/coin"
                    fi
                fi
                # PKG_LIBS="-L/usr/local/lib -Wl,-rpath /usr/local/lib"
            else
                AC_MSG_ERROR([none of the directories /usr/local/lib, /usr/local/lib64, /usr/lib of /usr/lib64 exists])
            fi

        else
    
            AC_MSG_NOTICE([found executable clp in $prefix])

            if test [ -d "$prefix/include" ] ; then
                PKG_CPPFLAGS="-I$prefix/include"
                if test [ -d "$prefix/include/coin" ] ; then
                    PKG_CPPFLAGS="$PKG_CPPFLAGS -I$prefix/include/coin"
                fi
            else
                AC_MSG_ERROR([directory $prefix/include does not exist])
            fi

            if test [ -d "$prefix/lib" -o -d "$prefix/lib64" ] ; then
                PKG_LIBS=""
                if test [ -d "$prefix/lib" ] ; then
                    PKG_LIBS="$PKG_LIBS -L$prefix/lib"
                    if test [ -d "$prefix/lib/coin" ] ; then
                        PKG_LIBS="$PKG_LIBS -L$prefix/lib/coin"
                    fi
                fi
                if test [ -d "$prefix/lib64" ] ; then
                    PKG_LIBS="$PKG_LIBS -L$prefix/lib64"
                    if test [ -d "$prefix/lib64/coin" ] ; then
                        PKG_LIBS="$PKG_LIBS -L$prefix/lib64/coin"
                    fi
                fi
            else
                AC_MSG_ERROR([neither $prefix/lib nor $prefix/lib64 exist])
            fi

        fi
    
    fi

    PKG_LIBS="${PKG_LIBS} -lClp -lCoinUtils -lm"

    ARGS_OK="TRUE"
fi


dnl # --------------------------------------------------------------------
dnl # CFLAGS
dnl # --------------------------------------------------------------------

if test [ -z "$PKG_CFLAGS" ] ; then
    if test [ "${R_ARCH}" = "/x86_64" ] ; then
        PKG_CFLAGS="-m64"
    elif test [ "${R_ARCH}" = "/i386" ] ; then
        PKG_CFLAGS="-m32"
    else
        PKG_CFLAGS=""
    fi
fi


dnl # --------------------------------------------------------------------
dnl # check header and library
dnl # --------------------------------------------------------------------

LIBS="${PKG_LIBS}"
CFLAGS="${CFLAGS} ${PKG_CFLAGS}"
CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"

AC_CHECK_HEADERS([stdbool.h],,
    AC_MSG_ERROR([Could not find stdbool.h]))

AC_CHECK_HEADER([Clp_C_Interface.h],,
    AC_MSG_ERROR([Could not find Clp_C_Interface.h:
      pcaL1 requires clp from http://www.coin-or.org/projects/Clp.xml
      use --with-clp-include or CLP_INCLUDE to specify the include path.]),
      [#ifdef HAVE_STDBOOL_H
      # include <stdbool.h>
      #endif
      ])

AC_SEARCH_LIBS([Clp_newModel], [Clp], ,
    AC_MSG_ERROR([Could not link to clp:
      use --with-clp-lib or CLP_LIB to specify the lib path.]))



dnl # --------------------------------------------------------------------
dnl # check for new functions in Clp-1.17.2
dnl # --------------------------------------------------------------------

AC_SEARCH_LIBS([Clp_setRowName], [Clp], AC_DEFINE([HAVE_CLP_EXT1_17_2], [1], [defined, if extended functions are available]), AC_MSG_NOTICE([Clp extensions 1.17.2 not available]), )





dnl # --------------------------------------------------------------------
dnl # substitute src/Makevars
dnl # --------------------------------------------------------------------

AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(PKG_CPPFLAGS)
AC_OUTPUT(src/Makevars)

exit 0
