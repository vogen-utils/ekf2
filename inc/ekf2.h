#ifndef ekf2_h
#define ekf2_h

/** parameters */

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

/** paremeters end */

#include <ae2f/cc/qualifier.h>

#undef	ekf2_real_t
#if	EKF2_REAL_BIT <= 32
#	define	ekf2_real_t	float
#elif	EKF2_REAL_BIT <= 64
#	define	ekf2_real_t	double
#endif

#undef	ekf2_real_sfx
#if	EKF2_REAL_BIT <= 32
#	define	ekf2_real_sfx(a)	a##f
#elif	EKF2_REAL_BIT <= 64
#	define	ekf2_real_sfx(a)	a
#endif

enum EKF2_RET_ {
	EKF2_RET_OK,
	EKF2_RET_NIL_HANDLE,
	EKF2_RET_INIT_DIM_TOO_BIG,
	EKF2_RET_UPDATE_FAILED_INV_SYM
};

enum EKF2_DIM_MAX_ {
	EKF2_DIM_MAX_STATE		= 1 << EKF2_DIM_MAX_STATE_EXP2,
	EKF2_DIM_MAX_MEASUREMENT	= 1 << EKF2_DIM_MAX_MEASUREMENT_EXP2,
	EKF2_DIM_MAX_CTL		= 1 << EKF2_DIM_MAX_CTL_EXP2
};

typedef ekf2_real_t ekf2_vsreal_t __attribute__ ((vector_size (sizeof(ekf2_real_t) * (EKF2_DIM_MAX_STATE))));
typedef ekf2_real_t ekf2_vmreal_t __attribute__ ((vector_size (sizeof(ekf2_real_t) * (EKF2_DIM_MAX_MEASUREMENT))));
typedef ekf2_real_t ekf2_vcreal_t __attribute__ ((vector_size (sizeof(ekf2_real_t) * (EKF2_DIM_MAX_CTL))));

typedef	struct {
	ekf2_vsreal_t(*m_func)(
			const ekf2_vsreal_t* const ae2f_restrict	rd_x,
			const ekf2_vcreal_t* const ae2f_restrict	rd_u_opt,
			const ekf2_real_t				c_dt,
			ae2f_unused void*				h_usr
		     );

	void(*m_jac)(
			const ekf2_vsreal_t* const ae2f_restrict	rd_x,
			const ekf2_vcreal_t* const ae2f_restrict	rd_u_opt,
			const ekf2_real_t				c_dt,
			ekf2_vsreal_t (* const ae2f_restrict		ret)[EKF2_DIM_MAX_STATE],
			ae2f_unused void*				h_usr
		    );
} ekf2_func_state_transition_t;

typedef	struct {
	ekf2_vmreal_t(*m_func)(
			const ekf2_vsreal_t* const ae2f_restrict	rd_x,
			void*					h_usr
		     );

	void(*m_jac)(
			const ekf2_vsreal_t* const ae2f_restrict	rd_x,
			ekf2_vsreal_t (* const ae2f_restrict	r_z)[EKF2_DIM_MAX_MEASUREMENT],
			void*					h_usr
		    );
} ekf2_func_measurement_t;

typedef struct {
	void*				m_usr_handle;
	ekf2_func_state_transition_t	m_func_state_transition;
	ekf2_func_measurement_t		m_func_measurement;

	ekf2_vsreal_t	m_state_estimate;
	ekf2_vsreal_t	m_err_cov[EKF2_DIM_MAX_STATE];
	ekf2_vsreal_t	m_pnoise[EKF2_DIM_MAX_STATE];
	ekf2_vmreal_t	m_mnoise[EKF2_DIM_MAX_MEASUREMENT];

	ekf2_dim_t	m_dim_state;
	ekf2_dim_t	m_dim_measurement;
	ekf2_dim_t	m_dim_ctl;
} ekf2_ctx_t;

enum EKF2_RET_ init_ekf2(
		ekf2_ctx_t* ae2f_restrict h_ekf,
		ekf2_dim_t	c_dim_state,
		ekf2_dim_t	c_dim_measure,
		ekf2_dim_t	c_dim_ctl,
		ekf2_func_state_transition_t fn_state_transition,
		ekf2_func_measurement_t fn_measurement
		);

enum EKF2_RET_ ekf2_predict(ekf2_ctx_t* ae2f_restrict h_ekf2, const ekf2_vcreal_t* ae2f_restrict const rd_u, const ekf2_real_t c_dt);
enum EKF2_RET_ ekf2_update(ekf2_ctx_t* ae2f_restrict h_ekf2, const ekf2_vmreal_t *ae2f_restrict const rd_z);

#endif
