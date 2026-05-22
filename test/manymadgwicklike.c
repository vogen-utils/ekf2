/*
 * Quaternion-based attitude EKF with gyro-bias estimation.
 *
 *   state x = [ q0  q1  q2  q3  bx  by  bz ]^T    (7)
 *   control u = [ wx  wy  wz ]^T                  (3, raw gyro)
 *   measurement z = [ ax  ay  az ]^T              (3, accelerometer)
 *
 * Dynamics:  q_dot = 0.5 * q (x) (0, w - b),    b_dot = 0
 * Measurement: h(x) = R(q)^T * [0,0,1]   (gravity in body frame, +Z up)
 *
 * Gyro feeds in as the control input; accelerometer is the measurement.
 * EKF2_DIM_MAX_CTL_EXP2 is bumped to 2 so the control vector can hold 3 elements.
 */

#include <ekf2.h>
#include <stdio.h>
#include <math.h>

enum { IDX_Q0 = 0, IDX_Q1, IDX_Q2, IDX_Q3, IDX_BX, IDX_BY, IDX_BZ };

#define DIM_STATE 7
#define DIM_MEAS  3
#define DIM_CTL   3

/* ---- state transition f(x, u, dt) ---- */
static ekf2_vsreal_t f_state(
		const ekf2_vsreal_t* const ae2f_restrict	rd_x,
		const ekf2_vcreal_t* const ae2f_restrict	rd_u,
		const ekf2_real_t				dt,
		ae2f_unused void*				h_usr)
{
	const ekf2_real_t q0 = (*rd_x)[IDX_Q0];
	const ekf2_real_t q1 = (*rd_x)[IDX_Q1];
	const ekf2_real_t q2 = (*rd_x)[IDX_Q2];
	const ekf2_real_t q3 = (*rd_x)[IDX_Q3];
	const ekf2_real_t bx = (*rd_x)[IDX_BX];
	const ekf2_real_t by = (*rd_x)[IDX_BY];
	const ekf2_real_t bz = (*rd_x)[IDX_BZ];

	const ekf2_real_t wx = (*rd_u)[0] - bx;
	const ekf2_real_t wy = (*rd_u)[1] - by;
	const ekf2_real_t wz = (*rd_u)[2] - bz;

	const ekf2_real_t h = ekf2_real_sfx(0.5) * dt;

	ekf2_vsreal_t out = *rd_x;
	out[IDX_Q0] += h * (-q1*wx - q2*wy - q3*wz);
	out[IDX_Q1] += h * ( q0*wx + q2*wz - q3*wy);
	out[IDX_Q2] += h * ( q0*wy - q1*wz + q3*wx);
	out[IDX_Q3] += h * ( q0*wz + q1*wy - q2*wx);

	/* renormalise to keep |q| = 1 */
	const ekf2_real_t n = ekf2_real_sfx(sqrt)(
			  out[IDX_Q0]*out[IDX_Q0] + out[IDX_Q1]*out[IDX_Q1]
			+ out[IDX_Q2]*out[IDX_Q2] + out[IDX_Q3]*out[IDX_Q3]);
	const ekf2_real_t inv = ekf2_real_sfx(1.) / n;
	out[IDX_Q0] *= inv;
	out[IDX_Q1] *= inv;
	out[IDX_Q2] *= inv;
	out[IDX_Q3] *= inv;

	/* bias rows are carried straight through (already in `out` from copy of rd_x) */
	return out;
}

/* ---- state Jacobian F = df/dx, evaluated at the prior x ---- */
static void f_jac(
		const ekf2_vsreal_t* const ae2f_restrict	rd_x,
		const ekf2_vcreal_t* const ae2f_restrict	rd_u,
		const ekf2_real_t				dt,
		ekf2_vsreal_t (* const ae2f_restrict		ret)[EKF2_DIM_MAX_STATE],
		ae2f_unused void*				h_usr)
{
	const ekf2_real_t q0 = (*rd_x)[IDX_Q0];
	const ekf2_real_t q1 = (*rd_x)[IDX_Q1];
	const ekf2_real_t q2 = (*rd_x)[IDX_Q2];
	const ekf2_real_t q3 = (*rd_x)[IDX_Q3];
	const ekf2_real_t bx = (*rd_x)[IDX_BX];
	const ekf2_real_t by = (*rd_x)[IDX_BY];
	const ekf2_real_t bz = (*rd_x)[IDX_BZ];

	const ekf2_real_t wx = (*rd_u)[0] - bx;
	const ekf2_real_t wy = (*rd_u)[1] - by;
	const ekf2_real_t wz = (*rd_u)[2] - bz;

	const ekf2_real_t h = ekf2_real_sfx(0.5) * dt;
	const ekf2_real_t one = ekf2_real_sfx(1.);

	/* caller zeros `ret`, so only non-zero entries need to be written */

	/* dq/dq:  I + (dt/2) * Omega(wc) */
	(*ret)[0][0] = one;   (*ret)[0][1] = -h*wx; (*ret)[0][2] = -h*wy; (*ret)[0][3] = -h*wz;
	(*ret)[1][0] =  h*wx; (*ret)[1][1] = one;   (*ret)[1][2] =  h*wz; (*ret)[1][3] = -h*wy;
	(*ret)[2][0] =  h*wy; (*ret)[2][1] = -h*wz; (*ret)[2][2] = one;   (*ret)[2][3] =  h*wx;
	(*ret)[3][0] =  h*wz; (*ret)[3][1] =  h*wy; (*ret)[3][2] = -h*wx; (*ret)[3][3] = one;

	/* dq/db = -dq/dwc  (since wc = u - b) */
	(*ret)[0][IDX_BX] =  h*q1; (*ret)[0][IDX_BY] =  h*q2; (*ret)[0][IDX_BZ] =  h*q3;
	(*ret)[1][IDX_BX] = -h*q0; (*ret)[1][IDX_BY] =  h*q3; (*ret)[1][IDX_BZ] = -h*q2;
	(*ret)[2][IDX_BX] = -h*q3; (*ret)[2][IDX_BY] = -h*q0; (*ret)[2][IDX_BZ] =  h*q1;
	(*ret)[3][IDX_BX] =  h*q2; (*ret)[3][IDX_BY] = -h*q1; (*ret)[3][IDX_BZ] = -h*q0;

	/* db/db = I,  db/dq = 0 */
	(*ret)[IDX_BX][IDX_BX] = one;
	(*ret)[IDX_BY][IDX_BY] = one;
	(*ret)[IDX_BZ][IDX_BZ] = one;
}

/* ---- measurement h(x) = R(q)^T * [0,0,1] ---- */
static ekf2_vmreal_t h_meas(
		const ekf2_vsreal_t* const ae2f_restrict	rd_x,
		ae2f_unused void*				h_usr)
{
	const ekf2_real_t q0 = (*rd_x)[IDX_Q0];
	const ekf2_real_t q1 = (*rd_x)[IDX_Q1];
	const ekf2_real_t q2 = (*rd_x)[IDX_Q2];
	const ekf2_real_t q3 = (*rd_x)[IDX_Q3];

	ekf2_vmreal_t out = (ekf2_vmreal_t){0,};
	out[0] = ekf2_real_sfx(2.) * (q1*q3 + q0*q2);
	out[1] = ekf2_real_sfx(2.) * (q2*q3 - q0*q1);
	out[2] = ekf2_real_sfx(1.) - ekf2_real_sfx(2.) * (q1*q1 + q2*q2);
	return out;
}

/* ---- measurement Jacobian H = dh/dx ---- */
static void h_jac(
		const ekf2_vsreal_t* const ae2f_restrict	rd_x,
		ekf2_vsreal_t (* const ae2f_restrict		ret)[EKF2_DIM_MAX_MEASUREMENT],
		ae2f_unused void*				h_usr)
{
	const ekf2_real_t q0 = (*rd_x)[IDX_Q0];
	const ekf2_real_t q1 = (*rd_x)[IDX_Q1];
	const ekf2_real_t q2 = (*rd_x)[IDX_Q2];
	const ekf2_real_t q3 = (*rd_x)[IDX_Q3];

	const ekf2_real_t two  = ekf2_real_sfx(2.);
	const ekf2_real_t four = ekf2_real_sfx(4.);

	/* caller zeros `ret`; bias columns stay 0 (h doesn't depend on bias) */

	(*ret)[0][IDX_Q0] =  two*q2;
	(*ret)[0][IDX_Q1] =  two*q3;
	(*ret)[0][IDX_Q2] =  two*q0;
	(*ret)[0][IDX_Q3] =  two*q1;

	(*ret)[1][IDX_Q0] = -two*q1;
	(*ret)[1][IDX_Q1] = -two*q0;
	(*ret)[1][IDX_Q2] =  two*q3;
	(*ret)[1][IDX_Q3] =  two*q2;

	(*ret)[2][IDX_Q0] =  ekf2_real_sfx(0.);
	(*ret)[2][IDX_Q1] = -four*q1;
	(*ret)[2][IDX_Q2] = -four*q2;
	(*ret)[2][IDX_Q3] =  ekf2_real_sfx(0.);
}

int main(void);
int main(void) {
	ekf2_ctx_t ctx = {0};

	const ekf2_func_state_transition_t f_st = { .m_func = f_state, .m_jac = f_jac };
	const ekf2_func_measurement_t      f_z  = { .m_func = h_meas,  .m_jac = h_jac  };

	if(init_ekf2(&ctx, DIM_STATE, DIM_MEAS, DIM_CTL, f_st, f_z) != EKF2_RET_OK)
		return 1;

	/* initial state: identity quaternion, zero bias */
	ctx.m_state_estimate[IDX_Q0] = ekf2_real_sfx(1.);

	/* initial covariance P (diagonal) */
	ctx.m_err_cov[IDX_Q0][IDX_Q0] = ekf2_real_sfx(1e-2);
	ctx.m_err_cov[IDX_Q1][IDX_Q1] = ekf2_real_sfx(1e-2);
	ctx.m_err_cov[IDX_Q2][IDX_Q2] = ekf2_real_sfx(1e-2);
	ctx.m_err_cov[IDX_Q3][IDX_Q3] = ekf2_real_sfx(1e-2);
	ctx.m_err_cov[IDX_BX][IDX_BX] = ekf2_real_sfx(1e-3);
	ctx.m_err_cov[IDX_BY][IDX_BY] = ekf2_real_sfx(1e-3);
	ctx.m_err_cov[IDX_BZ][IDX_BZ] = ekf2_real_sfx(1e-3);

	/* process noise Q: quaternion gets gyro-driven noise, bias is a slow random walk */
	ctx.m_pnoise[IDX_Q0][IDX_Q0] = ekf2_real_sfx(1e-6);
	ctx.m_pnoise[IDX_Q1][IDX_Q1] = ekf2_real_sfx(1e-6);
	ctx.m_pnoise[IDX_Q2][IDX_Q2] = ekf2_real_sfx(1e-6);
	ctx.m_pnoise[IDX_Q3][IDX_Q3] = ekf2_real_sfx(1e-6);
	ctx.m_pnoise[IDX_BX][IDX_BX] = ekf2_real_sfx(1e-8);
	ctx.m_pnoise[IDX_BY][IDX_BY] = ekf2_real_sfx(1e-8);
	ctx.m_pnoise[IDX_BZ][IDX_BZ] = ekf2_real_sfx(1e-8);

	/* accelerometer noise R */
	ctx.m_mnoise[0][0] = ekf2_real_sfx(1e-2);
	ctx.m_mnoise[1][1] = ekf2_real_sfx(1e-2);
	ctx.m_mnoise[2][2] = ekf2_real_sfx(1e-2);

	/* simulate stationary, level: gyro = 0, accel reads unit gravity along +Z */
	const ekf2_real_t dt = ekf2_real_sfx(0.01);
	ekf2_vcreal_t gyro  = (ekf2_vcreal_t){0,};
	ekf2_vmreal_t accel = (ekf2_vmreal_t){ sin(1), 0, cos(1) };

	for(int k = 0; k < 2000000000; ++k) {
		ekf2_predict(&ctx, &gyro, dt);
		ekf2_update(&ctx, &accel);
	}

	printf("final quaternion: [% .4f % .4f % .4f % .4f]\n",
			(double)ctx.m_state_estimate[IDX_Q0],
			(double)ctx.m_state_estimate[IDX_Q1],
			(double)ctx.m_state_estimate[IDX_Q2],
			(double)ctx.m_state_estimate[IDX_Q3]);
	printf("final gyro bias:  [% .4f % .4f % .4f]\n",
			(double)ctx.m_state_estimate[IDX_BX],
			(double)ctx.m_state_estimate[IDX_BY],
			(double)ctx.m_state_estimate[IDX_BZ]);

	return 0;
}
