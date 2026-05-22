#include "ekf2.h"
#include "./helpers.h"
#include <ae2f/cc/branches/naive.h>

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

	for(ekf2_dim_t i = EKF2_DIM_MAX_STATE; i --> 0; ) {
		h_ekf->m_pnoise[i] = h_ekf->m_err_cov[i] = (ekf2_vsreal_t){0, };
	}

	for(ekf2_dim_t i = EKF2_DIM_MAX_MEASUREMENT; i --> 0; ) {
		h_ekf->m_mnoise[i] = (ekf2_vmreal_t) {0, };
	}

	return EKF2_RET_OK;
}

enum EKF2_RET_ ekf2_predict(ekf2_ctx_t* ae2f_restrict h_ekf, const ekf2_vcreal_t* ae2f_restrict const rd_u, const ekf2_real_t c_dt)
{
	ekf2_vsreal_t	STATE_TRANS_JACRET[EKF2_DIM_MAX_STATE];

	ae2f_ifezerr(h_ekf)	return EKF2_RET_NIL_HANDLE;

	for(ekf2_dim_t C = EKF2_DIM_MAX_STATE; C --> 0; ) {
		STATE_TRANS_JACRET[C] = (ekf2_vsreal_t) { 0, };
	}

	h_ekf->m_func_state_transition.m_jac(
			&h_ekf->m_state_estimate
			, rd_u, c_dt, &STATE_TRANS_JACRET
			, h_ekf->m_usr_handle);

	h_ekf->m_state_estimate = h_ekf->m_func_state_transition.m_func(
			&h_ekf->m_state_estimate
			, rd_u, c_dt
			, h_ekf->m_usr_handle
			);

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

enum EKF2_RET_ ekf2_update(ekf2_ctx_t* ae2f_restrict h_ekf, const ekf2_vmreal_t* ae2f_restrict const rd_z)
{
	ekf2_vmreal_t	Y;
	ekf2_vsreal_t	H[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	S[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	S_INV[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	PHt[EKF2_DIM_MAX_STATE] = {0, };
	ekf2_vmreal_t	K[EKF2_DIM_MAX_STATE] = {0, };
	ekf2_vsreal_t	IKH[EKF2_DIM_MAX_STATE] = {0, };
	ekf2_vsreal_t	P_NEW[EKF2_DIM_MAX_STATE] = {0, };

	ae2f_ifezerr(h_ekf)	return EKF2_RET_NIL_HANDLE;
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

	ae2f_ifnezerr_strict(s_mat_inv_sym(&S, &S_INV, h_ekf->m_dim_measurement))
		return EKF2_RET_UPDATE_FAILED_INV_SYM;

	s_mat_PHT(
			&h_ekf->m_err_cov
			, &H
			, &PHt
			, h_ekf->m_dim_measurement
			, h_ekf->m_dim_state
		 );

	for(ekf2_dim_t i = h_ekf->m_dim_state; i --> 0; ) {
		for(ekf2_dim_t k = h_ekf->m_dim_measurement; k --> 0; )
			for(ekf2_dim_t j = h_ekf->m_dim_measurement; j --> 0; )
				K[i][j] += PHt[i][k] * S_INV[k][j];
	}


	for(ekf2_dim_t i = h_ekf->m_dim_state; i --> 0; ) {
		const ekf2_vmreal_t Ky = K[i] * Y;
		for(ekf2_dim_t j = h_ekf->m_dim_measurement; j --> 0; )
			h_ekf->m_state_estimate[i] += Ky[j];
	}

	for(ekf2_dim_t i = h_ekf->m_dim_state; i --> 0; ) IKH[i][i] = ekf2_real_sfx(1.);
	for(ekf2_dim_t i = h_ekf->m_dim_state; i --> 0; ) {
		for(ekf2_dim_t k = h_ekf->m_dim_measurement; k --> 0; )
			for(ekf2_dim_t j = h_ekf->m_dim_state; j --> 0; )
				IKH[i][j] -= K[i][k] * H[k][j];
	}

	s_matmul_ss(
			&IKH
			, &h_ekf->m_err_cov
			, &P_NEW
			, h_ekf->m_dim_state
		   );

	for(ekf2_dim_t i = EKF2_DIM_MAX_STATE; i --> h_ekf->m_dim_state; )
		h_ekf->m_err_cov[i] = P_NEW[i];

	for(ekf2_dim_t i = 0; i < h_ekf->m_dim_state; ++i) {
		h_ekf->m_err_cov[i][i] = P_NEW[i][i];
		for(ekf2_dim_t j = i + 1; j < h_ekf->m_dim_state; ++j) {
			const ekf2_real_t avg = (P_NEW[i][j] + P_NEW[j][i])
				* ekf2_real_sfx(0.5);
			h_ekf->m_err_cov[i][j] = avg;
			h_ekf->m_err_cov[j][i] = avg;
		}
	}

	return EKF2_RET_OK;
}
