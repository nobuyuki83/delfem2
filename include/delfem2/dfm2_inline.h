#ifdef DFM2_INLINE
#undef DFM2_INLINE
#endif

#ifdef DFM2_HEADER_ONLY
#  define DFM2_INLINE inline
#else
#  define DFM2_INLINE
#endif
