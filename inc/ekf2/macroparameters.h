#ifndef	ekf2_dim_t
typedef	unsigned	ekf2_dim_t;
#	define	ekf2_dim_t	ekf2_dim_t
#endif

#ifndef	ekf2_bool_t
typedef	unsigned	ekf2_bool_t;
#	define ekf2_bool_t	ekf2_bool_t
#endif

#ifndef	EKF2_REAL_BIT
#	define	EKF2_REAL_BIT	32
#endif

#ifndef	EKF2_DIM_MAX_STATE_EXP2
#	define	EKF2_DIM_MAX_STATE_EXP2	3
#endif

#ifndef	EKF2_DIM_MAX_MEASUREMENT_EXP2
#	define	EKF2_DIM_MAX_MEASUREMENT_EXP2	2
#endif

#ifndef	EKF2_DIM_MAX_CTL_EXP2
#	define	EKF2_DIM_MAX_CTL_EXP2	2
#endif

#ifndef	EKF2_DIM_MAX_STATE
#	define	EKF2_DIM_MAX_STATE		(1 << EKF2_DIM_MAX_STATE_EXP2)
#endif

#ifndef	EKF2_DIM_MAX_CTL
#	define	EKF2_DIM_MAX_CTL		(1 << EKF2_DIM_MAX_CTL_EXP2)
#endif

#ifndef	EKF2_DIM_MAX_MEASUREMENT
#	define	EKF2_DIM_MAX_MEASUREMENT	(1 << EKF2_DIM_MAX_MEASUREMENT_EXP2)
#endif

#ifndef	EKF2_HEADER_ONLY
#	define	EKF2_HEADER_ONLY	0
#endif
