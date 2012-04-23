AC_INIT([pcaL1], [1.1], [jpbrooks@vcu.edu])

dnl # --------------------------------------------------------------------
dnl # R compiler and flags
dnl # --------------------------------------------------------------------

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "could not determine R_HOME"
    exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC


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


AC_ARG_ENABLE(lapack,
              AC_HELP_STRING([--enable-lapack[=yes|no]],
                             [link to liblapack and libblas [[yes]]]
                            ),
              [LAPACK=$enableval], [LAPACK="yes"])


dnl # --------------------------------------------------------------------
dnl # test arguments
dnl # --------------------------------------------------------------------

INCL="-lClp -lCoinUtils -lm -lbz2 -lz"
if test "${LAPACK}" != "no"; then
    INCL="$INCL -llapack -lblas"
fi


dnl # check environment variables for include and library directories
if test [ -n "$CLP_INCLUDE" -a -n "$CLP_LIB" ] ; then
     AC_MSG_NOTICE([trying directories from environment variables])
    if test [ -d "$CLP_INCLUDE" ] ; then
        CLP_CPPFLAGS="-I${CLP_INCLUDE}/coin"
    else
        AC_MSG_ERROR([directory $CLP_INCLUDE does not exist])
    fi

    if test [ -d "$CLP_LIB" ] ; then
        if test [ -d "$CLP_LIB/coin" ] ; then
            LIBS="-L${CLP_LIB} -L${CLP_LIB}/coin $INCL"
        else
            LIBS="-L${CLP_LIB} $INCL"
        fi
    else
        AC_MSG_ERROR([directory $CLP_LIB does not exist])
    fi

else

dnl # if no arguments are given, try pkg-config.  If that doesn't work,
dnl # try to find the clp executable. If is not
dnl # found, use /usr/local


    if test [ "$prefix" = "NONE" ] ; then

        AC_MSG_NOTICE([No prefix given, trying pkg-config])
	AC_PATH_PROG([PKGCONF],[pkg-config],[],[$PATH:/usr/local/bin:ext/bin:ext:/sw/bin:/opt/bin])
	if test "xx$PKGCONF" != "xx"; then
	    AC_MSG_CHECKING([whether pkg-config knows about Clp])
	    if "${PKGCONF}" --exists clp; then
		AC_MSG_RESULT([yes])
		LIBS=`${PKGCONF} --libs clp`
		CLP_CPPFLAGS=`${PKGCONF} --cflags-only-I clp`
	    else 
		AC_MSG_NOTICE([pkg-config did not find it, trying to find path to clp])
		AC_PREFIX_PROGRAM(clp)
		if test [ "$prefix" = "NONE" ] ; then
		    AC_MSG_NOTICE([path not found, trying -I/usr/local/include and -L/usr/local/lib])
        	    if test [ -d "/usr/local/include" -a -d "/usr/local/lib" ] ; then
        	        LIBS="-L/usr/local/lib -L/usr/local/lib/coin $INCL"
	     	        CLP_CPPFLAGS="-I/usr/local/include/coin"
        	    else
        	        AC_MSG_ERROR([directories /usr/local/include and /usr/local/lib do not exist])
        	    fi
		else
		    if test [ -d "$prefix/include" -a -d "$prefix/lib" ] ; then
        	        CLP_CPPFLAGS="-I$prefix/include -I$prefix/include/coin"
        	        LIBS="-L$prefix/lib -L$prefix/lib/coin $INCL"
        	    else
        	        AC_MSG_ERROR([directories $prefix/include and $prefix/lib do not exist])
        	    fi
		fi
	    fi
	else
	    AC_MSG_NOTICE([pkg-config not found, trying to find path to clp])
	    AC_PREFIX_PROGRAM(clp)
	    if test [ "$prefix" = "NONE" ] ; then
	        AC_MSG_NOTICE([path not found, trying -I/usr/local/include and -L/usr/local/lib])
                if test [ -d "/usr/local/include" -a -d "/usr/local/lib" ] ; then
                    LIBS="-L/usr/local/lib -L/usr/local/lib/coin $INCL"
	            CLP_CPPFLAGS="-I/usr/local/include -I/usr/local/include/coin"
                else
                    AC_MSG_ERROR([directories /usr/local/include and /usr/local/lib do not exist])
                fi
	    fi
	fi
    else
	AC_MSG_NOTICE([trying prefix])

        if test [ -d "$prefix/include" -a -d "$prefix/lib" ] ; then
            CLP_CPPFLAGS="-I$prefix/include -I$prefix/include/coin"
            LIBS="-L$prefix/lib -L$prefix/lib/coin $INCL"
        else
            AC_MSG_ERROR([directories $prefix/include and $prefix/lib do not exist])
        fi

    fi

fi


CPPFLAGS="${CPPFLAGS} ${CLP_CPPFLAGS}"

dnl # --------------------------------------------------------------------
dnl # check header and library
dnl # --------------------------------------------------------------------

AC_CHECK_HEADER([Clp_C_Interface.h],,
    AC_MSG_ERROR([Could not find Clp_C_Interface.h:
      pcaL1 requires clp from http://www.coin-or.org/projects/Clp.xml
      use --with-clp-include or CLP_INCLUDE to specify the include path.]))


AC_SEARCH_LIBS([Clp], ,
    AC_MSG_ERROR([Could not link to clp:
      use --with-clp-lib or CLP_LIB to specify the lib path.]))


dnl # --------------------------------------------------------------------
dnl # substitute src/Makevars
dnl # --------------------------------------------------------------------

AC_SUBST(CLP_CPPFLAGS)
AC_SUBST(LIBS)
AC_OUTPUT(src/Makevars)

exit 0
