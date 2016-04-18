/*
 * Adapted from Open Shading Language with this license:
 *
 * Copyright (c) 2009-2010 Sony Pictures Imageworks Inc., et al.
 * All Rights Reserved.
 *
 * Modifications Copyright 2011, Blender Foundation.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * * Neither the name of Sony Pictures Imageworks nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __BSDF_DISNEY_DIFFUSE_
#define __BSDF_DISNEY_DIFFUSE_H__

CCL_NAMESPACE_BEGIN

/* DIFFUSE */

ccl_device int bsdf_disney_diffuse_setup(ShaderClosure *sc)
{
	sc->type = CLOSURE_BSDF_DISNEY_DIFFUSE_ID;
	return SD_BSDF|SD_BSDF_HAS_EVAL;
}

ccl_device float3 bsdf_disney_diffuse_eval_reflect(const ShaderClosure *sc,
        float3 color, const float3 I, const float3 omega_in, float *pdf)
{
	float3 N = sc->N;

	float cos_pi = fmaxf(dot(N, omega_in), 0.0f) * M_1_PI_F;
	*pdf = cos_pi;
	return make_float3(cos_pi * color[0], cos_pi * color[1], cos_pi * color[2]);
}

ccl_device float3 bsdf_disney_diffuse_eval_transmit(const ShaderClosure *sc,
        float3 color, const float3 I, const float3 omega_in, float *pdf)
{
	return make_float3(0.0f, 0.0f, 0.0f);
}

ccl_device int bsdf_disney_diffuse_sample(const ShaderClosure *sc, float3 color,
        float3 Ng, float3 I, float3 dIdx, float3 dIdy, float randu, float randv,
        float3 *eval, float3 *omega_in, float3 *domega_in_dx,
        float3 *domega_in_dy, float *pdf)
{
	float3 N = sc->N;

	// distribution over the hemisphere
	sample_cos_hemisphere(N, randu, randv, omega_in, pdf);

	if(dot(Ng, *omega_in) > 0.0f) {
		*eval = make_float3(*pdf * color[0], *pdf * color[1], *pdf * color[2]);
#ifdef __RAY_DIFFERENTIALS__
		// TODO: find a better approximation for the diffuse bounce
		*domega_in_dx = (2 * dot(N, dIdx)) * N - dIdx;
		*domega_in_dy = (2 * dot(N, dIdy)) * N - dIdy;
#endif
	}
	else
		*pdf = 0.0f;

	return LABEL_REFLECT|LABEL_DIFFUSE;
}

CCL_NAMESPACE_END

#endif /* __BSDF_DISNEY_DIFFUSE_H__ */
