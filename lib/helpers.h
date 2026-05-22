#ifndef	ekf2_helpers_h
#define	ekf2_helpers_h

#include <ekf2.h>
#include <ae2f/cc/inline.h>
#include <ae2f/cc/branches/strict.h>
#include <math.h>

extern ekf2_real_t	ekf2_real_sfx(sqrt)(const ekf2_real_t);

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
		for(ekf2_dim_t	j = c_dim_state; j --> 0; ) {
			const ekf2_vsreal_t	AB = (*rd_a)[i] * (*rd_b)[j];
			for(ekf2_dim_t k = c_dim_state; k --> 0; )
				T[i][j] += AB[k];
		}

	for(ekf2_dim_t	i = EKF2_DIM_MAX_STATE; i --> 0; )
		(*ret)[i] = T[i];
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
		for(ekf2_dim_t k = c_dim_state; k --> 0; ) {
			for(ekf2_dim_t j = c_dim_measure; j --> 0; )
				(*ret)[i][j] += (HP)[i][k] * (*rd_h)[j][k];
		}
}

ae2f_inline_f ekf2_bool_t	s_mat_inv_sym(
		const ekf2_vmreal_t (* const ae2f_restrict	rd_a)[EKF2_DIM_MAX_MEASUREMENT],
		ekf2_vmreal_t (* const ae2f_restrict		ret_a_inv)[EKF2_DIM_MAX_MEASUREMENT],
		const ekf2_dim_t				c_dim_measure
		) {
	ekf2_vmreal_t	L[EKF2_DIM_MAX_MEASUREMENT] = {0, };
	ekf2_vmreal_t	Y[EKF2_DIM_MAX_MEASUREMENT] = {0, };

	for(ekf2_dim_t i = 0; i < c_dim_measure; ++i) {
		for(ekf2_dim_t j = 0; j <= i; ++j) {
			ekf2_vmreal_t	SM	= L[i] * L[j];
			ekf2_real_t	S	= (*rd_a)[i][j];

			for(ekf2_dim_t k = j; k --> 0; )
				S -= SM[k];

			if(i ^ j) {
				L[i][j] = S / L[j][j];
			} else {
				ae2f_ifnezerr_strict(S <= ekf2_real_sfx(0.))	return 1;
				L[i][j]	= ekf2_real_sfx(sqrt)(S);
			}
		}	
	}

	for(ekf2_dim_t col = c_dim_measure; col --> 0; ) {
		for(ekf2_dim_t i = 0; i < c_dim_measure; ++i) {
			ekf2_real_t S = (ekf2_real_t)(i == col);
			for(ekf2_dim_t k = i; k --> 0; )
				S -= L[i][k] * Y[k][col];

			Y[i][col] = S / L[i][i];
		}
	}

	for(ekf2_dim_t i = EKF2_DIM_MAX_MEASUREMENT; i --> 0; )
		(*ret_a_inv)[i] = (ekf2_vmreal_t) {0, };


	for(ekf2_dim_t col = c_dim_measure; col --> 0; )
		for(ekf2_dim_t i = c_dim_measure; i-- > 0;) {
			ekf2_real_t S = Y[i][col];
			for(ekf2_dim_t k = i + 1; k < c_dim_measure; ++k) {
				S -= L[k][i] * (*ret_a_inv)[k][col];
			}
			(*ret_a_inv)[i][col] = S / L[i][i];
		}

	return 0;
}

ae2f_inline_f 	static void s_mat_PHT(
		const ekf2_vsreal_t (* ae2f_restrict const rd_p)[EKF2_DIM_MAX_STATE],
		const ekf2_vsreal_t (* ae2f_restrict const rd_h)[EKF2_DIM_MAX_MEASUREMENT],
		ekf2_vmreal_t (* ae2f_restrict const ret_k)[EKF2_DIM_MAX_STATE],
		const ekf2_dim_t		c_dim_measure,
		const ekf2_dim_t		c_dim_state
		) {
	for(ekf2_dim_t i = c_dim_state; i --> 0;) {
		(*ret_k)[i] = (ekf2_vmreal_t) {0, };
		for(ekf2_dim_t k = c_dim_state; k --> 0; ) {
			for(ekf2_dim_t j = c_dim_measure; j --> 0; )
				(*ret_k)[i][j] += (*rd_p)[i][k] * (*rd_h)[j][k];
		}
	}
}

#endif
