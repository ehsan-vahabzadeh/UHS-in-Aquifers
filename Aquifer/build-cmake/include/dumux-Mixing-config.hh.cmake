
#ifndef DUMUX_MIXING_CONFIG_HH
#define DUMUX_MIXING_CONFIG_HH

/* Define to 1 if you have module dumux-Mixing available */
#ifndef HAVE_DUMUX_MIXING
#cmakedefine01 HAVE_DUMUX_MIXING
#endif

/* begin dumux-hyuspre
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/



/* Define to the version of dumux-hyuspre */
#define DUMUX_HYUSPRE_VERSION "@DUMUX_HYUSPRE_VERSION@"

/* Define to the major version of dumux-hyuspre */
#define DUMUX_HYUSPRE_VERSION_MAJOR @DUMUX_HYUSPRE_VERSION_MAJOR@

/* Define to the minor version of dumux-hyuspre */
#define DUMUX_HYUSPRE_VERSION_MINOR @DUMUX_HYUSPRE_VERSION_MINOR@

/* Define to the revision of dumux-hyuspre */
#define DUMUX_HYUSPRE_VERSION_REVISION @DUMUX_HYUSPRE_VERSION_REVISION@

/* end dumux-hyuspre
   Everything below here will be overwritten
*/



#if __has_include(<dumux-config.hh>)
  #include <dumux-config.hh>
#endif


#endif // DUMUX_MIXING_CONFIG_HH


#ifndef DUMUX_MIXING_CONFIG_BOTTOM_HH
#define DUMUX_MIXING_CONFIG_BOTTOM_HH



#endif // DUMUX_MIXING_CONFIG_BOTTOM_HH
