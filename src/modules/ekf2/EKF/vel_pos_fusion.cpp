/****************************************************************************
 *
 *   Copyright (c) 2015 Estimation and Control Library (ECL). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name ECL nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file vel_pos_fusion.cpp
 *
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Siddharth Bharat Purohit <siddharthbharatpurohit@gmail.com>
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include <mathlib/mathlib.h>
#include "ekf.h"

void Ekf::updateVelocity(const uint64_t& sample_time_us, const Vector3f& velocity, const Vector3f& obs_var, const float innov_gate, estimator_aid_source_3d_s& vel_aid_src)
{
	resetEstimatorAidStatus(vel_aid_src);

	for (int i = 0; i < 3; i++) {
		vel_aid_src.observation[i] = velocity(i);
		vel_aid_src.observation_variance[i] = obs_var(i);

		vel_aid_src.innovation[i] = _state.vel(i) - velocity(i);
		vel_aid_src.innovation_variance[i] = P(4 + i, 4 + i) + obs_var(i);
	}

	setEstimatorAidStatusTestRatio(vel_aid_src, innov_gate);

	// vz special case if there is bad vertical acceleration data, then don't reject measurement,
	// but limit innovation to prevent spikes that could destabilise the filter
	if (_fault_status.flags.bad_acc_vertical && vel_aid_src.innovation_rejected[2]) {
		const float innov_limit = innov_gate * sqrtf(vel_aid_src.innovation_variance[2]);
		vel_aid_src.innovation[2] = math::constrain(vel_aid_src.innovation[2], -innov_limit, innov_limit);
		vel_aid_src.innovation_rejected[2] = false;
	}

	vel_aid_src.timestamp_sample = sample_time_us;
}

void Ekf::updatePosition(const uint64_t& sample_time_us, const Vector3f& position, const Vector3f& obs_var, const float innov_gate, estimator_aid_source_3d_s& pos_aid_src)
{
	resetEstimatorAidStatus(pos_aid_src);

	for (int i = 0; i < 3; i++) {
		pos_aid_src.observation[i] = position(i);
		pos_aid_src.observation_variance[i] = obs_var(i);

		pos_aid_src.innovation[i] = _state.pos(i) - position(i);
		pos_aid_src.innovation_variance[i] = P(7 + i, 7 + i) + obs_var(i);
	}

	setEstimatorAidStatusTestRatio(pos_aid_src, innov_gate);

	// z special case if there is bad vertical acceleration data, then don't reject measurement,
	// but limit innovation to prevent spikes that could destabilise the filter
	if (_fault_status.flags.bad_acc_vertical && pos_aid_src.innovation_rejected[2]) {
		const float innov_limit = innov_gate * sqrtf(pos_aid_src.innovation_variance[2]);
		pos_aid_src.innovation[2] = math::constrain(pos_aid_src.innovation[2], -innov_limit, innov_limit);
		pos_aid_src.innovation_rejected[2] = false;
	}

	pos_aid_src.timestamp_sample = sample_time_us;
}

void Ekf::fuseVelocity(estimator_aid_source_3d_s& vel_aid_src)
{
	// vx & vy
	if (vel_aid_src.fusion_enabled[0] && vel_aid_src.fusion_enabled[1] &&
		!vel_aid_src.innovation_rejected[0] && !vel_aid_src.innovation_rejected[1]) {
		for (int i = 0; i < 2; i++) {
			if (fuseVelPosHeight(vel_aid_src.innovation[i], vel_aid_src.innovation_variance[i], i)) {
				vel_aid_src.fused[i] = true;
				vel_aid_src.time_last_fuse[i] = _time_last_imu;
			}
		}
	}

	// vz
	if (vel_aid_src.fusion_enabled[2] && !vel_aid_src.innovation_rejected[2]) {
		if (fuseVelPosHeight(vel_aid_src.innovation[2], vel_aid_src.innovation_variance[2], 2)) {
			vel_aid_src.fused[2] = true;
			vel_aid_src.time_last_fuse[2] = _time_last_imu;
		}
	}
}

void Ekf::fusePosition(estimator_aid_source_3d_s& pos_aid_src)
{
	// x & y
	if (pos_aid_src.fusion_enabled[0] && !pos_aid_src.innovation_rejected[0]
	 && pos_aid_src.fusion_enabled[1] && !pos_aid_src.innovation_rejected[1]
	) {
		for (int i = 0; i < 2; i++) {
			if (fuseVelPosHeight(pos_aid_src.innovation[i], pos_aid_src.innovation_variance[i], 3 + i)) {
				pos_aid_src.fused[i] = true;
				pos_aid_src.time_last_fuse[i] = _time_last_imu;
			}
		}
	}

	// z
	if (pos_aid_src.fusion_enabled[2] && !pos_aid_src.innovation_rejected[2]) {
		if (fuseVelPosHeight(pos_aid_src.innovation[2], pos_aid_src.innovation_variance[2], 5)) {
			pos_aid_src.fused[2] = true;
			pos_aid_src.time_last_fuse[2] = _time_last_imu;
		}
	}
}

// Helper function that fuses a single velocity or position measurement
bool Ekf::fuseVelPosHeight(const float innov, const float innov_var, const int obs_index)
{
	Vector24f Kfusion;  // Kalman gain vector for any single observation - sequential fusion is used.
	const unsigned state_index = obs_index + 4;  // we start with vx and this is the 4. state

	// calculate kalman gain K = PHS, where S = 1/innovation variance
	for (int row = 0; row < _k_num_states; row++) {
		Kfusion(row) = P(row, state_index) / innov_var;
	}

	SquareMatrix24f KHP;

	for (unsigned row = 0; row < _k_num_states; row++) {
		for (unsigned column = 0; column < _k_num_states; column++) {
			KHP(row, column) = Kfusion(row) * P(state_index, column);
		}
	}

	// if the covariance correction will result in a negative variance, then
	// the covariance matrix is unhealthy and must be corrected
	bool healthy = true;

	for (int i = 0; i < _k_num_states; i++) {
		if (P(i, i) < KHP(i, i)) {
			// zero rows and columns
			P.uncorrelateCovarianceSetVariance<1>(i, 0.0f);

			healthy = false;
		}
	}

	setVelPosStatus(obs_index, healthy);

	if (healthy) {
		// apply the covariance corrections
		P -= KHP;

		fixCovarianceErrors(true);

		// apply the state corrections
		fuse(Kfusion, innov);

		return true;
	}

	return false;
}

void Ekf::setVelPosStatus(const int index, const bool healthy)
{
	switch (index) {
	case 0:
		if (healthy) {
			_fault_status.flags.bad_vel_N = false;
			_time_last_hor_vel_fuse = _time_last_imu;

		} else {
			_fault_status.flags.bad_vel_N = true;
		}

		break;

	case 1:
		if (healthy) {
			_fault_status.flags.bad_vel_E = false;
			_time_last_hor_vel_fuse = _time_last_imu;

		} else {
			_fault_status.flags.bad_vel_E = true;
		}

		break;

	case 2:
		if (healthy) {
			_fault_status.flags.bad_vel_D = false;
			_time_last_ver_vel_fuse = _time_last_imu;

		} else {
			_fault_status.flags.bad_vel_D = true;
		}

		break;

	case 3:
		if (healthy) {
			_fault_status.flags.bad_pos_N = false;
			_time_last_hor_pos_fuse = _time_last_imu;

		} else {
			_fault_status.flags.bad_pos_N = true;
		}

		break;

	case 4:
		if (healthy) {
			_fault_status.flags.bad_pos_E = false;
			_time_last_hor_pos_fuse = _time_last_imu;

		} else {
			_fault_status.flags.bad_pos_E = true;
		}

		break;

	case 5:
		if (healthy) {
			_fault_status.flags.bad_pos_D = false;
			_time_last_hgt_fuse = _time_last_imu;

		} else {
			_fault_status.flags.bad_pos_D = true;
		}

		break;
	}
}
