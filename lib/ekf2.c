#include "./ekf2.h"
#include <ae2f/cc/branches/naive.h>
#include <ae2f/cc/inline.h>

enum EKF2_RET_ init_ekf2(
		ekf2_ctx_t* ae2f_restrict h_ekf,
		ekf2_dim_t	c_dim_state,
		ekf2_dim_t	c_dim_measure,
		ekf2_dim_t	c_dim_ctl,
		ekf2_func_state_transition_t fn_state_transition,
		ekf2_func_measurement_t fn_measurement
		)
{
	ae2f_ifezerr(h_ekf)
		return EKF2_RET_NIL_HANDLE;
	ae2f_ifnezerr_strict(c_dim_state > EKF2_DIM_MAX_STATE)
		return EKF2_RET_INIT_DIM_TOO_BIG;
	ae2f_ifnezerr_strict(c_dim_measure > EKF2_DIM_MAX_MEASUREMENT)
		return EKF2_RET_INIT_DIM_TOO_BIG;
	ae2f_ifnezerr_strict(c_dim_ctl > EKF2_DIM_MAX_CTL)
		return EKF2_RET_INIT_DIM_TOO_BIG;

	h_ekf->m_dim_state		= c_dim_state;
	h_ekf->m_dim_measurement	= c_dim_measure;
	h_ekf->m_dim_ctl		= c_dim_ctl;

	h_ekf->m_func_state_transition	= fn_state_transition;
	h_ekf->m_func_measurement	= fn_measurement;

	h_ekf->m_state_estimate = (ekf2_vsreal_t){0, };

	return EKF2_RET_OK;
}

ae2f_inline_f	static void s_matmul_ss(
		ekf2_vsreal_t (* ae2f_restrict const		ret)[EKF2_DIM_MAX_STATE],
		const ekf2_vsreal_t (* ae2f_restrict const	rd_a)[EKF2_DIM_MAX_STATE],
		const ekf2_vsreal_t (* ae2f_restrict const	rd_b)[EKF2_DIM_MAX_STATE],
		const ekf2_dim_t				c_dim_state
		) {
	ekf2_vsreal_t	T[EKF2_DIM_MAX_STATE] = { 0, };

	for(ekf2_dim_t	i = c_dim_state; i --> 0; )
		for(ekf2_dim_t	j = c_dim_state; j --> 0; )
			for(ekf2_dim_t k = c_dim_state; k --> 0; )
				T[i][j] += (*rd_a)[i][k] * (*rd_b)[k][j];

	for(ekf2_dim_t	i = EKF2_DIM_MAX_STATE; i --> 0; )
		(*ret)[i] = T[i];
}

ae2f_inline_f	static void s_matmul_trans_ss(
		ekf2_vsreal_t (* ae2f_restrict const		ret)[EKF2_DIM_MAX_STATE],
		const ekf2_vsreal_t (* ae2f_restrict const	rd_a)[EKF2_DIM_MAX_STATE],
		const ekf2_vsreal_t (* ae2f_restrict const	rd_b)[EKF2_DIM_MAX_STATE],
		const ekf2_dim_t				c_dim_state
		)
{
	ekf2_vsreal_t	T[EKF2_DIM_MAX_STATE] = { 0, };

	for(ekf2_dim_t	i = c_dim_state; i --> 0; )
		for(ekf2_dim_t	j = c_dim_state; j --> 0; )
			for(ekf2_dim_t k = c_dim_state; k --> 0; )
				T[i][j] += (*rd_a)[i][k] * (*rd_b)[j][k];

	for(ekf2_dim_t	i = EKF2_DIM_MAX_STATE; i --> 0; )
		(*ret)[i] = T[i];
}

enum EKF2_RET_ ekf2_predict(ekf2_ctx_t* ae2f_restrict h_ekf, const ekf2_real_t* ae2f_restrict const rd_u, const ekf2_real_t c_dt)
{
	ekf2_vsreal_t	STATE_TRANS_JACRET[EKF2_DIM_MAX_STATE];

	ae2f_ifnezerr(h_ekf)	return EKF2_RET_NIL_HANDLE;

	h_ekf->m_state_estimate = h_ekf->m_func_state_transition.m_func(
			&h_ekf->m_state_estimate
			, rd_u, c_dt
			, h_ekf->m_usr_handle
			);

	for(ekf2_dim_t C = EKF2_DIM_MAX_STATE; C --> 0; ) {
		STATE_TRANS_JACRET[C] = (ekf2_vsreal_t) { 0, };
	}

	h_ekf->m_func_state_transition.m_jac(
			&h_ekf->m_state_estimate
			, rd_u, c_dt, &STATE_TRANS_JACRET
			, h_ekf->m_usr_handle);

	{
		ekf2_vsreal_t	JACRET_NOISE[EKF2_DIM_MAX_STATE];
		ekf2_vsreal_t	JACRET_NOISE_JACRET[EKF2_DIM_MAX_STATE];

		s_matmul_ss(
				&JACRET_NOISE
				, &STATE_TRANS_JACRET
				, &h_ekf->m_err_cov
				, h_ekf->m_dim_state);

		s_matmul_trans_ss(
				&JACRET_NOISE_JACRET
				, &JACRET_NOISE
				, &STATE_TRANS_JACRET
				, h_ekf->m_dim_state
				);

		for(ekf2_dim_t i = h_ekf->m_dim_state; i --> 0; )
			h_ekf->m_err_cov[i] = h_ekf->m_pnoise[i] + JACRET_NOISE_JACRET[i];
	}

	return EKF2_RET_OK;
}

ae2f_inline_f void s_mat_HPHT(
		const ekf2_vsreal_t(* ae2f_restrict const rd_h)[EKF2_DIM_MAX_MEASUREMENT],
		const ekf2_vsreal_t(* ae2f_restrict const rd_p)[EKF2_DIM_MAX_STATE],
		ekf2_vmreal_t(* ae2f_restrict const ret)[EKF2_DIM_MAX_MEASUREMENT],
		const ekf2_dim_t		c_dim_measure,
		const ekf2_dim_t		c_dim_state
		)
{
	ekf2_vsreal_t HP[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	for(ekf2_dim_t	i = c_dim_measure; i --> 0; )
		for(ekf2_dim_t k = c_dim_state; k --> 0; )
			for(ekf2_dim_t j = c_dim_state; j --> 0; )
				HP[i][j] += (*rd_h)[i][k] * (*rd_p)[k][j];

	for(ekf2_dim_t i = EKF2_DIM_MAX_MEASUREMENT; i --> 0; )
		(*ret)[i] = (ekf2_vmreal_t) { 0, };

	for(ekf2_dim_t	i = c_dim_measure; i --> 0; )
		for(ekf2_dim_t k = c_dim_state; k --> 0; )
			for(ekf2_dim_t j = c_dim_measure; j --> 0; )
				(*ret)[i][j] += (HP)[i][k] * (*rd_h)[j][k];
}

extern ekf2_real_t	ekf2_real_sfx(sqrt)(const ekf2_real_t);

ae2f_inline_f ekf2_bool_t	s_mat_inv_sym(
		const ekf2_vmreal_t (* const ae2f_restrict	rd_a)[EKF2_DIM_MAX_MEASUREMENT],
		const ekf2_vmreal_t (* const ae2f_restrict	ret_a_inv)[EKF2_DIM_MAX_MEASUREMENT],
		const ekf2_dim_t				c_dim_measure
		) {
	ekf2_vmreal_t	L[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	Y[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	for(ekf2_dim_t i = c_dim_measure; i --> 0; )
		for(ekf2_dim_t j = c_dim_measure; j --> 0; ) {
			ekf2_vmreal_t	SM	= L[i] * L[j];
			ekf2_real_t	S	= (*rd_a)[i][j];

			for(ekf2_dim_t k = j; j --> 0; )
				S -= SM[k];

			if(i ^ j) {
				L[i][j] = S / L[i][j];
			} else {
				if(S <= ekf2_real_sfx(0.))	return 1;
				L[i][j]	= ekf2_real_sfx(sqrt)(S);
			}
		}

	for(ekf2_dim_t col = c_dim_measure; col --> 0; ) {
		ekf2_vmreal_t SV = 
	}
}


enum EKF2_RET_ ekf2_update(ekf2_ctx_t* ae2f_restrict h_ekf, const ekf2_vmreal_t* ae2f_restrict const rd_z)
{
	ekf2_vmreal_t	Y;
	ekf2_vsreal_t	H[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	S[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	S_INV[EKF2_DIM_MAX_MEASUREMENT] = {0, };

	ae2f_ifnezerr(h_ekf)	return EKF2_RET_NIL_HANDLE;
	Y = *rd_z - h_ekf->m_func_measurement.m_func(
			&h_ekf->m_state_estimate
			, h_ekf->m_usr_handle
			);

	h_ekf->m_func_measurement.m_jac(&h_ekf->m_state_estimate, &H, h_ekf->m_usr_handle);
	s_mat_HPHT(
			&H
			, &h_ekf->m_err_cov
			, &S
			, h_ekf->m_dim_measurement
			, h_ekf->m_dim_state
			);

	for(ekf2_dim_t i = h_ekf->m_dim_measurement; i --> 0; )
		S[i] += h_ekf->m_mnoise[i];

	
}
