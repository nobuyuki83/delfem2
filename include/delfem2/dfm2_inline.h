#ifdef DFM2_INLINE
#  undef DFM2_INLINE
#endif

#ifndef DFM2_STATIC_LIBRARY
#  define DFM2_INLINE inline
#else
#  define DFM2_INLINE
#endif
