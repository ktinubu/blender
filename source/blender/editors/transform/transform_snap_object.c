/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/transform/transform_snap_object.c
 *  \ingroup edtransform
 */

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_kdopbvh.h"
#include "BLI_memarena.h"
#include "BLI_ghash.h"
#include "BLI_linklist.h"
#include "BLI_listbase.h"
#include "BLI_utildefines.h"

#include "DNA_armature_types.h"
#include "DNA_curve_types.h"
#include "DNA_scene_types.h"
#include "DNA_object_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_screen_types.h"
#include "DNA_view3d_types.h"

#include "BKE_DerivedMesh.h"
#include "BKE_object.h"
#include "BKE_anim.h"  /* for duplis */
#include "BKE_editmesh.h"
#include "BKE_main.h"
#include "BKE_tracking.h"

#include "ED_transform.h"
#include "ED_transform_snap_object_context.h"
#include "ED_view3d.h"
#include "ED_armature.h"

#include "transform.h"

enum eViewProj {
	VIEW_PROJ_NONE = -1,
	VIEW_PROJ_ORTHO = 0,
	VIEW_PROJ_PERSP = -1,
};

/* Flags related to occlusion planes */
enum {
	BEHIND_A_PLANE      = 0,
	ISECT_CLIP_PLANE    = 1 << 0,
	IN_FRONT_ALL_PLANES = 1 << 1,
	TEST_RANGE_DEPTH    = 1 << 2,
};

typedef struct SnapData {
	float ray_origin[3];
	float ray_start[3];
	float ray_dir[3];
	float pmat[4][4]; /* perspective matrix */

	short mval[2];
	float win_half[2];/* win x and y */
	float depth_range[2];/* zNear and zFar of Perspective View */

	short snap_to_flag;
	enum eViewProj view_proj;
	bool test_occlusion;

	float (*rv3d_clip)[4];
	float (*dinamic_plane)[4];

} SnapData;

typedef struct SnapPorjectVert {
	short co[2];
	float depth;
	bool is_visible;
} SnapPorjectVert;

typedef struct SnapObjectData {
	enum {
		SNAP_MESH = 1,
		SNAP_EDIT_MESH,
	} type;

	SnapPorjectVert *proj_vert;

} SnapObjectData;

typedef struct SnapDataMesh {
	SnapObjectData sd;

	BVHTreeFromMesh treedata;
	/* extend */
	MPoly *mpoly;
	bool poly_allocated;
	bool has_looptris;

	float imat[4][4];

} SnapDataMesh;

typedef struct SnapDataEditMesh {
	SnapObjectData sd;

	BVHTreeFromEditMesh treedata;
	BMEditMesh *em;
	float imat[4][4];

} SnapDataEditMesh;

/** used for storing objects in scene */
typedef struct SnapObjectBase {
	struct SnapObjectBase *next, *prev;

	struct Object *ob;
	bool is_obedit;
	float mat[4][4];
	void *data;
} SnapObjectBase;

struct SnapObjectContext {
	Main *bmain;
	Scene *scene;
	int flag;

	/* Optional: when performing screen-space projection.
	 * otherwise this doesn't take viewport into account. */
	bool use_v3d;
	struct {
		const struct View3D *v3d;
		const struct ARegion *ar;
	} v3d_data;

	/* Object -> SnapObjectData map */
	struct {
		ListBase object_map;
		MemArena *mem_arena;
	} cache;

	/* Filter data, returns true to check this value */
	struct {
		struct {
			bool (*test_vert_fn)(BMVert *, void *user_data);
			bool (*test_edge_fn)(BMEdge *, void *user_data);
			bool (*test_face_fn)(BMFace *, void *user_data);
			void *user_data;
		} edit_mesh;
	} callbacks;

};


/* -------------------------------------------------------------------- */

/** \name Support for storing all depths, not just the first (raycast 'all')
 *
 * This uses a list of #SnapObjectHitDepth structs.
 *
 * \{ */

/* Store all ray-hits */
struct RayCastAll_Data {
	void *bvhdata;

	/* internal vars for adding depths */
	BVHTree_RayCastCallback raycast_callback;

	float(*obmat)[4];
	float(*timat)[3];

	float len_diff;
	float local_scale;

	Object *ob;
	unsigned int ob_uuid;

	/* output data */
	ListBase *hit_list;
	bool retval;
};

static struct SnapObjectHitDepth *hit_depth_create(
        const float depth, const float co[3], const float no[3], int index,
        Object *ob, const float obmat[4][4], unsigned int ob_uuid)
{
	struct SnapObjectHitDepth *hit = MEM_mallocN(sizeof(*hit), __func__);

	hit->depth = depth;
	copy_v3_v3(hit->co, co);
	copy_v3_v3(hit->no, no);
	hit->index = index;

	hit->ob = ob;
	copy_m4_m4(hit->obmat, (float(*)[4])obmat);
	hit->ob_uuid = ob_uuid;

	return hit;
}

static int hit_depth_cmp_cb(const void *arg1, const void *arg2)
{
	const struct SnapObjectHitDepth *h1 = arg1;
	const struct SnapObjectHitDepth *h2 = arg2;
	int val = 0;

	if (h1->depth < h2->depth) {
		val = -1;
	}
	else if (h1->depth > h2->depth) {
		val = 1;
	}

	return val;
}

static void raycast_all_cb(void *userdata, int index, const BVHTreeRay *ray, BVHTreeRayHit *hit)
{
	struct RayCastAll_Data *data = userdata;
	data->raycast_callback(data->bvhdata, index, ray, hit);
	if (hit->index != -1) {
		/* get all values in worldspace */
		float location[3], normal[3];
		float depth;

		/* worldspace location */
		mul_v3_m4v3(location, data->obmat, hit->co);
		depth = (hit->dist + data->len_diff) / data->local_scale;

		/* worldspace normal */
		copy_v3_v3(normal, hit->no);
		mul_m3_v3(data->timat, normal);
		normalize_v3(normal);

		/* currently unused, and causes issues when looptri's haven't been calculated.
		 * since theres some overhead in ensuring this data is valid, it may need to be optional. */
#if 0
		if (data->dm) {
			hit->index = dm_looptri_to_poly_index(data->dm, &data->dm_looptri[hit->index]);
		}
#endif

		struct SnapObjectHitDepth *hit_item = hit_depth_create(
		        depth, location, normal, hit->index,
		        data->ob, data->obmat, data->ob_uuid);
		BLI_addtail(data->hit_list, hit_item);
	}
}

/** \} */


/* -------------------------------------------------------------------- */

/** \Common utilities
 * \{ */


MINLINE float depth_get(const float co[3], const float ray_start[3], const float ray_dir[3])
{
	float dvec[3];
	sub_v3_v3v3(dvec, co, ray_start);
	return dot_v3v3(dvec, ray_dir);
}

MINLINE void aabb_get_near_far_from_plane(
        const float plane_no[3],
        const float bbmin[3], const float bbmax[3],
        float bb_near[3], float bb_afar[3])
{
	if (plane_no[0] < 0) {
		bb_near[0] = bbmax[0];
		bb_afar[0] = bbmin[0];
	}
	else {
		bb_near[0] = bbmin[0];
		bb_afar[0] = bbmax[0];
	}
	if (plane_no[1] < 0) {
		bb_near[1] = bbmax[1];
		bb_afar[1] = bbmin[1];
	}
	else {
		bb_near[1] = bbmin[1];
		bb_afar[1] = bbmax[1];
	}
	if (plane_no[2] < 0) {
		bb_near[2] = bbmax[2];
		bb_afar[2] = bbmin[2];
	}
	else {
		bb_near[2] = bbmin[2];
		bb_afar[2] = bbmax[2];
	}
}

static bool snp_calc_project_vert(
        float (*clip_planes)[4], short plane_num, float pmat[4][4],
        const float win_half[2], bool is_persp, const float v[3],
        SnapPorjectVert *proj_vert)
{
	if (clip_planes) {
		if (!isect_point_planes_v3_negate(clip_planes, plane_num, v)) {
			// printf("snap_point behind clip plane\n");
			proj_vert->is_visible = false;
			return false;
		}
	}
	float pv[3] = {
	    (dot_m4_v3_row_x(pmat, v) + pmat[3][0]),
	    (dot_m4_v3_row_y(pmat, v) + pmat[3][1]),
	    (dot_m4_v3_row_z(pmat, v) + pmat[3][2]),
	};

	if (is_persp) {
		if (pv[2] < 0) {
			proj_vert->is_visible = false;
			return false;
		}
		mul_v3_fl(pv, 1 / mul_project_m4_v3_zfac(pmat, v));
	}

	if (fabs(pv[2] > 1)) {
		proj_vert->is_visible = false;
		return false;
	}

	proj_vert->co[0] = (short)(win_half[0] * (1.0f + pv[0]));
	proj_vert->co[1] = (short)(win_half[1] * (1.0f + pv[1]));
	proj_vert->depth = pv[2];
	proj_vert->is_visible = true;
	return true;
}

static float dist_aabb_to_plane(
        const float bbmin[3], const float bbmax[3],
        const float plane_co[3], const float plane_no[3])
{
	const float bb_near[3] = {
		(plane_no[0] < 0) ? bbmax[0] : bbmin[0],
		(plane_no[1] < 0) ? bbmax[1] : bbmin[1],
		(plane_no[2] < 0) ? bbmax[2] : bbmin[2],
	};
	return depth_get(bb_near, plane_co, plane_no);
}

/**
 * Check if a AABB is:
 * - BEHIND_A_PLANE     (0),
 * - ISECT_CLIP_PLANE   (1),
 * - IN_FRONT_ALL_PLANES(2)
*/
static short snp_isect_aabb_planes_v3(
        float(*planes)[4], int totplane, const float bbmin[3], const float bbmax[3])
{
	short ret = IN_FRONT_ALL_PLANES;

	float bb_near[3], bb_afar[3];
	for (int i = 0; i < totplane; i++) {
		aabb_get_near_far_from_plane(planes[i], bbmin, bbmax, bb_near, bb_afar);
		if (plane_point_side_v3(planes[i], bb_afar) < 0.0f) {
			return BEHIND_A_PLANE;
		}
		else if ((ret != ISECT_CLIP_PLANE) && (plane_point_side_v3(planes[i], bb_near) < 0.0f)) {
			ret = ISECT_CLIP_PLANE;
		}
	}

	return ret;
}

static void snp_nearest_local_data_get(
        SnapData *snpdt, float obmat[4][4],
        float (*local_pmat)[4], float (**local_clip_plane)[4], short *clip_plane_num)
{
	BLI_assert(local_clip_plane != NULL);

	mul_m4_m4m4(local_pmat, snpdt->pmat, obmat);;

	*clip_plane_num = 0;
	if (snpdt->dinamic_plane) {
		*clip_plane_num = 1;
	}
	if (snpdt->rv3d_clip) {
		*clip_plane_num += 4;
	}
	else if (!snpdt->dinamic_plane) {
		*local_clip_plane = NULL;
		return;
	}

	*local_clip_plane = MEM_mallocN(sizeof(**local_clip_plane) * (*clip_plane_num), __func__);

	float tobmat[4][4];
	transpose_m4_m4(tobmat, obmat);

	mul_v4_m4v4(*local_clip_plane[0], tobmat, snpdt->dinamic_plane[0]);

	for (short i = *clip_plane_num - 1; i != 0; i--) {
		mul_v4_m4v4(*local_clip_plane[i], tobmat, snpdt->rv3d_clip[i]);

//		mul_v4_fl(clip_local[i], 1 / len_v3(clip_local[i])); /* normalize plane */
	}
}

/**
 * Generates a struct with the immutable parameters that will be used on all objects.
 *
 * \param snap_to: Element to snap, Vertice, Edge or Face.
 * \param view_proj: ORTHO or PERSP.
 * Currently only works one at a time, but can eventually operate as flag.
 *
 * \param mval: Mouse coords.
 */
static bool snapdata_init_v3d(
        SnapData *snpdt, SnapObjectContext *sctx,
        const unsigned short snap_to_flag, const float mval[2], float *depth)
{
	if (!sctx->use_v3d) {
		return false;
	}

	snpdt->snap_to_flag = snap_to_flag;

	const ARegion *ar = sctx->v3d_data.ar;
	RegionView3D *rv3d = (RegionView3D *)ar->regiondata;

	snpdt->mval[0] = (short)mval[0];
	snpdt->mval[1] = (short)mval[1];

	ED_view3d_win_to_origin(ar, mval, snpdt->ray_origin);
	ED_view3d_win_to_vector(ar, mval, snpdt->ray_dir);

	ED_view3d_clip_range_get(
	        sctx->v3d_data.v3d, rv3d, &snpdt->depth_range[0], &snpdt->depth_range[1], false);

	madd_v3_v3v3fl(snpdt->ray_start, snpdt->ray_origin, snpdt->ray_dir, snpdt->depth_range[0]);

	if (rv3d->rflag & RV3D_CLIPPING) {
		snpdt->rv3d_clip = rv3d->clip;

		/* Get a new ray_start and a new depth to use in occlusion or snap to faces */
		float dummy_ray_end[3];
		madd_v3_v3v3fl(dummy_ray_end, snpdt->ray_origin, snpdt->ray_dir, snpdt->depth_range[1]);

		if (!clip_segment_v3_plane_n(
		        snpdt->ray_start, dummy_ray_end, snpdt->rv3d_clip, 4,
		        snpdt->ray_start, dummy_ray_end))
		{
			return false;
		}
		*depth = depth_get(dummy_ray_end, snpdt->ray_start, snpdt->ray_dir);
	}
	else {
		snpdt->rv3d_clip = NULL;
	}

	copy_m4_m4(snpdt->pmat, rv3d->persmat);
	snpdt->win_half[0] = ar->winx / 2;
	snpdt->win_half[1] = ar->winy / 2;

	snpdt->view_proj = rv3d->is_persp ? VIEW_PROJ_PERSP : VIEW_PROJ_ORTHO;
	snpdt->test_occlusion = true;

	snpdt->dinamic_plane = NULL;

	return true;
}

/**
 * Generates a struct with the immutable parameters that will be used on all objects.
 * Used only in ray_cast (snap to faces)
 * (ray-casting is handled without any projection matrix correction.)
 *
 * \param ray_origin: ray_start before being moved toward the ray_normal at the distance from vew3d clip_min.
 * \param ray_start: ray_origin moved for the start clipping plane (clip_min).
 * \param ray_direction: Unit length direction of the ray.
 * \param depth_range: distances of clipe plane min and clip plane max;
 */
static bool snapdata_init_ray(
        SnapData *snpdt,
        const float ray_start[3], const float ray_normal[3])
{
	snpdt->snap_to_flag = SCE_SELECT_FACE;

	copy_v3_v3(snpdt->ray_origin, ray_start);
	copy_v3_v3(snpdt->ray_start, ray_start);
	copy_v3_v3(snpdt->ray_dir, ray_normal);
//	snpdt->depth_range[0] = 0.0f;
//	snpdt->depth_range[1] = BVH_RAYCAST_DIST_MAX;

	snpdt->view_proj = VIEW_PROJ_NONE;
	snpdt->test_occlusion = false;

	return true;
}

static bool snap_point_v3(
        const short mval[2], const float co[3],
        float pmat[4][4], const float win_half[2], const bool is_persp,
        float(*planes)[4], const short totplane,
        short dist_px[2], float r_co[3])
{
	SnapPorjectVert proj_vert;
	if (!snp_calc_project_vert(
	        planes, totplane, pmat, win_half,
	        is_persp, co, &proj_vert))
	{
		return false;
	}

	short tmp_dist[2];
	tmp_dist[0] = abs(mval[0] - proj_vert.co[0]);
	tmp_dist[1] = abs(mval[1] - proj_vert.co[1]);

	if (tmp_dist[0] <= dist_px[0] && tmp_dist[1] <= dist_px[1]) {
		copy_v3_v3(r_co, co);
		copy_v2_v2_short(dist_px, tmp_dist);
		return true;
	}

	return false;
}

static bool snap_segment_v2v2(
        short snap_to, const short mval[2],
        SnapPorjectVert *proj_vert1, SnapPorjectVert *proj_vert2,
        short dist_2d[2], float *r_lambda)
{
	float r_close[2];
	short tmp_dist[2];
	float lambda = closest_to_line_segment_v2(
	        r_close, (float[2]){UNPACK2(mval)}, (float[2]){UNPACK2(proj_vert1->co)}, (float[2]){UNPACK2(proj_vert2->co)});

	if ((snap_to & SCE_SELECT_VERTEX) && (lambda < 0.25f || 0.75f < lambda)) {
		lambda = lambda > 0.5f;
		short *v = lambda ? proj_vert2->co : proj_vert1->co;
		tmp_dist[0] = abs(mval[0] - v[0]);
		tmp_dist[1] = abs(mval[1] - v[1]);

		if (tmp_dist[0] <= dist_2d[0] && tmp_dist[1] <= dist_2d[1]) {
			copy_v2_v2_short(dist_2d, tmp_dist);
			*r_lambda = lambda;
			return true;
		}
	}
	else {
		tmp_dist[0] = abs(mval[0] - (short)r_close[0]);
		tmp_dist[1] = abs(mval[1] - (short)r_close[1]);

		if (tmp_dist[0] <= dist_2d[0] && tmp_dist[1] <= dist_2d[1]) {
			copy_v2_v2_short(dist_2d, tmp_dist);
			*r_lambda = lambda * (proj_vert1->depth / proj_vert2->depth);
			return true;
		}
	}

	return false;
}

static bool snap_segment_v3v3(
        short snap_to, float (*clip_plane)[4], short clip_plane_num,
        float pmat[4][4], const short mval[2],
        const float win_half[2], const bool is_persp,
        const float va[3], const float vb[3],
        short dist_px[2], float r_co[3])
{
	bool ret = false;

	if (snap_to & SCE_SELECT_EDGE) {
		/* TODO (Germano): Compensate object scale */
		SnapPorjectVert proj_vert[2];
		snp_calc_project_vert(
		        clip_plane, clip_plane_num, pmat, win_half,
		        is_persp, va, &proj_vert[0]);
		snp_calc_project_vert(
		        clip_plane, clip_plane_num, pmat, win_half,
		        is_persp, vb, &proj_vert[1]);

		float lambda;

		if (snap_segment_v2v2(
		        snap_to, mval, &proj_vert[0], &proj_vert[1],
		        dist_px, &lambda))
		{
			interp_v3_v3v3(r_co, va, vb, lambda);
			return true;
		}
	}
	else {
		ret = snap_point_v3(
		        mval, va, pmat, win_half, is_persp,
		        clip_plane, clip_plane_num, dist_px, r_co);
		ret |= snap_point_v3(
		        mval, vb, pmat, win_half, is_persp,
		        clip_plane, clip_plane_num, dist_px, r_co);
	}

	return ret;
}

typedef struct SnapNearest2dPrecalc {
	struct {
		float ray_orig[3];
		float ray_dir[3];
		float pmat[4][4]; /* perspective matrix */
	} local;

	float ray_inv_dir[3];
	bool is_persp;
	float win_half[2];

	float mval[2];

	float depth_range[2];

} SnapNearest2dPrecalc;

static void snp_dist_squared_to_projected_aabb_precalc(
        struct SnapNearest2dPrecalc *nearest_precalc,
        float ray_orig[3], float ray_dir[3], float pmat[4][4],
        SnapData *snpdt)
{
	copy_v3_v3(nearest_precalc->local.ray_orig, ray_orig);
	copy_v3_v3(nearest_precalc->local.ray_dir, ray_dir);
	copy_m4_m4(nearest_precalc->local.pmat, pmat);
//	memcpy(&nearest_precalc->local, localdata, sizeof(nearest_precalc->local));

	for (int i = 0; i < 3; i++) {
		nearest_precalc->ray_inv_dir[i] =
		        (nearest_precalc->local.ray_dir[i] != 0.0f) ?
		        (1.0f / nearest_precalc->local.ray_dir[i]) : FLT_MAX;
	}

	nearest_precalc->is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
	copy_v2_v2(nearest_precalc->win_half, snpdt->win_half);
	copy_v2_v2(nearest_precalc->depth_range, snpdt->depth_range);

	copy_v2_v2(nearest_precalc->mval, (float[2]){UNPACK2(snpdt->mval)});
}

/* Returns the distance from a 2d coordinate to a BoundBox (Projected) */
static float snp_dist_squared_to_projected_aabb(
        struct SnapNearest2dPrecalc *data,
        const float bbmin[3], const float bbmax[3],
        short *flag)
{
	float bb_near[3], bb_afar[3];
	aabb_get_near_far_from_plane(
	        data->ray_inv_dir, bbmin, bbmax, bb_near, bb_afar);

	/* ISECT_CLIP_PLANE can replace this ? */
	if (*flag & TEST_RANGE_DEPTH) {
		/* Test if the entire AABB is behind us */
		float depth_near;
		float depth_afar;
		if (data->is_persp) {
			/* TODO: in perspective mode depth is not necessaly taken from bb_near and bb_afar */
			depth_near = mul_project_m4_v3_zfac(data->local.pmat, bb_near);
			depth_afar = mul_project_m4_v3_zfac(data->local.pmat, bb_afar);
			if (depth_afar < data->depth_range[0]) {
//				printf("Is behind near clip plane\n");
				return FLT_MAX;
			}
			if (depth_near > data->depth_range[1]) {
//				printf("Is after far clip plane\n");
				return FLT_MAX;
			}
			if (data->depth_range[0] < depth_near && depth_afar < data->depth_range[1]) {
				*flag &= ~TEST_RANGE_DEPTH;
			}
		}
		else {
			depth_near = dot_m4_v3_row_z(data->local.pmat, bb_near);
			depth_afar = dot_m4_v3_row_z(data->local.pmat, bb_afar);
			if (depth_afar < -1.0f) {
//				printf("Is behind near clip plane\n");
				return FLT_MAX;
			}
			if (depth_near > 1.0f) {
//				printf("Is after far clip plane\n");
				return FLT_MAX;
			}
			if (-1.0f < depth_near && depth_afar < 1.0f) {
				*flag &= ~TEST_RANGE_DEPTH;
			}
		}
	}

	const float tmin[3] = {
		(bb_near[0] - data->local.ray_orig[0]) * data->ray_inv_dir[0],
		(bb_near[1] - data->local.ray_orig[1]) * data->ray_inv_dir[1],
		(bb_near[2] - data->local.ray_orig[2]) * data->ray_inv_dir[2],
	};
	const float tmax[3] = {
		(bb_afar[0] - data->local.ray_orig[0]) * data->ray_inv_dir[0],
		(bb_afar[1] - data->local.ray_orig[1]) * data->ray_inv_dir[1],
		(bb_afar[2] - data->local.ray_orig[2]) * data->ray_inv_dir[2],
	};
	/* `va` and `vb` are the coordinates of the AABB edge closest to the ray */
	float va[3], vb[3];
	/* `rtmin` and `rtmax` are the minimum and maximum distances of the ray hits on the AABB */
	float rtmin, rtmax;
	int main_axis, tmin_axis, tmax_axis;

	if ((tmax[0] <= tmax[1]) && (tmax[0] <= tmax[2])) {
		rtmax = tmax[0];
		va[0] = vb[0] = bb_afar[0];
		tmin_axis = 0;
	}
	else if ((tmax[1] <= tmax[0]) && (tmax[1] <= tmax[2])) {
		rtmax = tmax[1];
		va[1] = vb[1] = bb_afar[1];
		tmin_axis = 1;
	}
	else {
		rtmax = tmax[2];
		va[2] = vb[2] = bb_afar[2];
		tmin_axis = 2;
	}

	if ((tmin[0] >= tmin[1]) && (tmin[0] >= tmin[2])) {
		rtmin = tmin[0];
		va[0] = vb[0] = bb_near[0];
		tmax_axis = 0;
	}
	else if ((tmin[1] >= tmin[0]) && (tmin[1] >= tmin[2])) {
		rtmin = tmin[1];
		va[1] = vb[1] = bb_near[1];
		tmax_axis = 1;
	}
	else {
		rtmin = tmin[2];
		va[2] = vb[2] = bb_near[2];
		tmax_axis = 2;
	}

	/* ray intersect `AABB` */
	if ((rtmin <= rtmax) || (tmin_axis == tmax_axis)) {
		return 0.0f;
	}

	main_axis = 3 - (tmin_axis + tmax_axis);

	va[main_axis] = bbmin[main_axis];
	vb[main_axis] = bbmax[main_axis];

	float scale = bbmax[main_axis] - bbmin[main_axis];

	float va2d[2] = {
		(dot_m4_v3_row_x(data->local.pmat, va) + data->local.pmat[3][0]),
		(dot_m4_v3_row_y(data->local.pmat, va) + data->local.pmat[3][1]),
	};
	float vb2d[2] = {
		(va2d[0] + data->local.pmat[main_axis][0] * scale),
		(va2d[1] + data->local.pmat[main_axis][1] * scale),
	};

	if (data->is_persp) {
		float depth_a = mul_project_m4_v3_zfac(data->local.pmat, va);
		float depth_b = depth_a + data->local.pmat[main_axis][3] * scale;
		va2d[0] /= depth_a;
		va2d[1] /= depth_a;
		vb2d[0] /= depth_b;
		vb2d[1] /= depth_b;
	}

	va2d[0] += 1.0f;
	va2d[1] += 1.0f;
	vb2d[0] += 1.0f;
	vb2d[1] += 1.0f;

	va2d[0] *= data->win_half[0];
	va2d[1] *= data->win_half[1];
	vb2d[0] *= data->win_half[0];
	vb2d[1] *= data->win_half[1];

	short dvec[2] = {data->mval[0] - va2d[0], data->mval[1] - va2d[1]};
	short edge[2] = {vb2d[0] - va2d[0], vb2d[1] - va2d[1]};
	float lambda = dvec[0] * edge[0] + dvec[1] * edge[1];
	if (lambda != 0.0f) {
		lambda /= edge[0] * edge[0] + edge[1] * edge[1];
		if (lambda <= 0.0f) {
			return len_squared_v2v2(data->mval, va2d);
		}
		else if (lambda >= 1.0f) {
			return len_squared_v2v2(data->mval, vb2d);
		}
		else {
			va2d[0] += edge[0] * lambda;
			va2d[1] += edge[1] * lambda;
			return len_squared_v2v2(data->mval, va2d);
		}
	}

	/* parallel ray */
	return len_squared_v2v2(data->mval, va2d);
}

static bool snap_boundbox_nearest_test(
        SnapData *snpdt, BoundBox *bb, float obmat[4][4], float dist_px)
{
	float bbmin[3], bbmax[3];
	mul_v3_m4v3(bbmin, obmat, bb->vec[0]);
	mul_v3_m4v3(bbmax, obmat, bb->vec[6]);

	if (snpdt->dinamic_plane) {
		if (snp_isect_aabb_planes_v3(
			        snpdt->dinamic_plane, 1,
			        bbmin, bbmax) == BEHIND_A_PLANE)
		{
			return false;
		}
	}
	if (snpdt->rv3d_clip) {
		if (snp_isect_aabb_planes_v3(
			        snpdt->rv3d_clip, 4,
			        bbmin, bbmax) == BEHIND_A_PLANE)
		{
			return false;
		}
	}

	struct SnapNearest2dPrecalc data;

	snp_dist_squared_to_projected_aabb_precalc(&data, snpdt->ray_origin, snpdt->ray_dir, snpdt->pmat, snpdt);

	short flag = TEST_RANGE_DEPTH;
	return snp_dist_squared_to_projected_aabb(
	        &data, bbmin, bbmax, &flag) < SQUARE(dist_px);
}

static bool snp_snap_boundbox_raycast_test(
        SnapData *snpdt, BoundBox *bb, float obmat[4][4])
{
	float bbmin[3], bbmax[3];
	mul_v3_m4v3(bbmin, obmat, bb->vec[0]);
	mul_v3_m4v3(bbmax, obmat, bb->vec[6]);

	/* was BKE_boundbox_ray_hit_check, see: cf6ca226fa58 */
	return isect_ray_aabb_v3_simple(
	        snpdt->ray_start, snpdt->ray_dir, bbmin, bbmax, NULL, NULL);
}

/** \} */


/* -------------------------------------------------------------------- */

/** \Utilities for DerivedMeshes and EditMeshes
* \{ */

static DerivedMesh *object_dm_final_get(Scene *scn, Object *ob)
{
	DerivedMesh *dm;

	/* in this case we want the mesh from the editmesh, avoids stale data. see: T45978 */
	BMEditMesh *em = BKE_editmesh_from_object(ob);
	if (em) {
		editbmesh_get_derived_cage_and_final(scn, ob, em, CD_MASK_BAREMESH, &dm);
	}
	else {
		dm = mesh_get_derived_final(scn, ob, CD_MASK_BAREMESH);
	}

	return dm;
}

/** \} */

/* -------------------------------------------------------------------- */

/** \name Internal Object Snapping API
 * \{ */

static bool snapArmature(
        SnapObjectContext *sctx, SnapData *snpdt,
        Object *ob, bArmature *arm, float obmat[4][4],
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float *UNUSED(r_no))
{
	bool retval = false;

	if (snpdt->snap_to_flag == SCE_SELECT_FACE) { /* Currently only edge and vert */
		return retval;
	}

	BoundBox *bb = BKE_object_boundbox_get(ob);
	if (bb && !snap_boundbox_nearest_test(snpdt, bb, obmat, *dist_px)) {
		return retval;
	}

	float local_pmat[4][4], (*local_clip_plane)[4];
	short clip_plane_num;
	snp_nearest_local_data_get(snpdt, obmat, local_pmat, &local_clip_plane, &clip_plane_num);

	bool is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
	short dist_2d[2] = {*dist_px, *dist_px};

	if (arm->edbo) {
		for (EditBone *eBone = arm->edbo->first; eBone; eBone = eBone->next) {
			if (eBone->layer & arm->layer) {
				/* skip hidden or moving (selected) bones */
				if ((eBone->flag & (BONE_HIDDEN_A | BONE_ROOTSEL | BONE_TIPSEL)) == 0) {
					retval |= snap_segment_v3v3(
					        snpdt->snap_to_flag, local_clip_plane, clip_plane_num,
					        local_pmat, snpdt->mval,
					        snpdt->win_half, is_persp,
					        eBone->head, eBone->tail,
					        dist_2d, r_loc);
				}
			}
		}
	}
	else if (ob->pose && ob->pose->chanbase.first) {
		for (bPoseChannel *pchan = ob->pose->chanbase.first; pchan; pchan = pchan->next) {
			Bone *bone = pchan->bone;
			/* skip hidden bones */
			if (bone && !(bone->flag & (BONE_HIDDEN_P | BONE_HIDDEN_PG))) {
				const float *head_vec = pchan->pose_head;
				const float *tail_vec = pchan->pose_tail;

				retval |= snap_segment_v3v3(
				        snpdt->snap_to_flag, local_clip_plane, clip_plane_num,
				        local_pmat, snpdt->mval,
				        snpdt->win_half, is_persp,
				        head_vec, tail_vec,
				        dist_2d, r_loc);
			}
		}
	}

	if (local_clip_plane) {
		MEM_freeN(local_clip_plane);
	}

	if (retval) {
		*dist_px = sqrt(dist_2d[0] * dist_2d[0] + dist_2d[1] * dist_2d[1]);
		mul_m4_v3(obmat, r_loc);
		*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);
		return true;
	}

	return false;
}

static bool snapCurve(
        SnapObjectContext *sctx, SnapData *snpdt,
        Object *ob, Curve *cu, float obmat[4][4],
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float *UNUSED(r_no))
{
	bool retval = false;

	/* only vertex snapping mode (eg control points and handles) supported for now) */
	if ((snpdt->snap_to_flag & SCE_SELECT_VERTEX) == 0) {
		return retval;
	}

	float local_pmat[4][4], (*local_clip_plane)[4];
	short clip_plane_num;
	snp_nearest_local_data_get(snpdt, obmat, local_pmat, &local_clip_plane, &clip_plane_num);
	bool is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
	short dist_2d[2] = {*dist_px, *dist_px};

	for (Nurb *nu = (ob->mode == OB_MODE_EDIT ? cu->editnurb->nurbs.first : cu->nurb.first); nu; nu = nu->next) {
		for (int u = 0; u < nu->pntsu; u++) {
			if (ob->mode == OB_MODE_EDIT) {
				if (nu->bezt) {
					/* don't snap to selected (moving) or hidden */
					if (nu->bezt[u].f2 & SELECT || nu->bezt[u].hide != 0) {
						break;
					}
					retval |= snap_point_v3(
					        snpdt->mval, nu->bezt[u].vec[1],
					        local_pmat, snpdt->win_half, is_persp,
					        local_clip_plane, clip_plane_num,
					        dist_2d, r_loc);
					/* don't snap if handle is selected (moving), or if it is aligning to a moving handle */
					if (!(nu->bezt[u].f1 & SELECT) &&
					    !(nu->bezt[u].h1 & HD_ALIGN && nu->bezt[u].f3 & SELECT))
					{
						retval |= snap_point_v3(
						        snpdt->mval, nu->bezt[u].vec[0],
						        local_pmat, snpdt->win_half, is_persp,
						        local_clip_plane, clip_plane_num,
						        dist_2d, r_loc);
					}
					if (!(nu->bezt[u].f3 & SELECT) &&
					    !(nu->bezt[u].h2 & HD_ALIGN && nu->bezt[u].f1 & SELECT))
					{
						retval |= snap_point_v3(
						        snpdt->mval, nu->bezt[u].vec[2],
						        local_pmat, snpdt->win_half, is_persp,
						        local_clip_plane, clip_plane_num,
						        dist_2d, r_loc);
					}
				}
				else {
					/* don't snap to selected (moving) or hidden */
					if (nu->bp[u].f1 & SELECT || nu->bp[u].hide != 0) {
						break;
					}
					retval |= snap_point_v3(
					        snpdt->mval, nu->bp[u].vec,
					        local_pmat, snpdt->win_half, is_persp,
					        local_clip_plane, clip_plane_num,
					        dist_2d, r_loc);
				}
			}
			else {
				/* curve is not visible outside editmode if nurb length less than two */
				if (nu->pntsu > 1) {
					if (nu->bezt) {
						retval |= snap_point_v3(
						        snpdt->mval, nu->bezt[u].vec[1],
						        local_pmat, snpdt->win_half, is_persp,
						        local_clip_plane, clip_plane_num,
						        dist_2d, r_loc);
					}
					else {
						retval |= snap_point_v3(
						        snpdt->mval, nu->bp[u].vec,
						        local_pmat, snpdt->win_half, is_persp,
						        local_clip_plane, clip_plane_num,
						        dist_2d, r_loc);
					}
				}
			}
		}
	}

	if (local_clip_plane) {
		MEM_freeN(local_clip_plane);
	}

	if (retval) {
		*dist_px = sqrt(dist_2d[0] * dist_2d[0] + dist_2d[1] * dist_2d[1]);
		mul_m4_v3(obmat, r_loc);
		*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);
		return true;
	}
	return false;
}

/* may extend later (for now just snaps to empty center) */
static bool snapEmpty(
        SnapObjectContext *sctx, SnapData *snpdt,
        Object *ob, float obmat[4][4],
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float *UNUSED(r_no))
{
	bool retval = false;

	if (ob->transflag & OB_DUPLI) {
		return retval;
	}

	/* for now only vertex supported */
	if (snpdt->snap_to_flag & SCE_SELECT_VERTEX) {
		bool is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
		short dist_2d[2] = {*dist_px, *dist_px};

		if (snpdt->dinamic_plane && !isect_point_planes_v3_negate(snpdt->dinamic_plane, 1, obmat[3])) {
			return false;
		}

		if (snap_point_v3(
			        snpdt->mval, obmat[3],
			        snpdt->pmat, snpdt->win_half, is_persp,
			        snpdt->rv3d_clip, 4,
			        dist_2d, r_loc))
		{
			*dist_px = sqrt(dist_2d[0] * dist_2d[0] + dist_2d[1] * dist_2d[1]);
			*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);
			retval = true;
		}
	}

	return retval;
}

static bool snapCamera(
        const SnapObjectContext *sctx, SnapData *snpdt,
        Object *object, float obmat[4][4],
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float *UNUSED(r_no))
{
	bool retval = false;

	if  (snpdt->snap_to_flag & SCE_SELECT_VERTEX) {
		Scene *scene = sctx->scene;

		MovieClip *clip = BKE_object_movieclip_get(scene, object, false);

		if (clip == NULL) {
			return retval;
		}
		if (object->transflag & OB_DUPLI) {
			return retval;
		}

		float orig_camera_mat[4][4], orig_camera_imat[4][4], imat[4][4];

		bool is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
		short dist_2d[2] = {*dist_px, *dist_px};

		BKE_tracking_get_camera_object_matrix(scene, object, orig_camera_mat);

		invert_m4_m4(orig_camera_imat, orig_camera_mat);
		invert_m4_m4(imat, obmat);

		RegionView3D *rv3d = sctx->v3d_data.ar->regiondata;

		MovieTracking *tracking = &clip->tracking;

		MovieTrackingObject *tracking_object;

		for (tracking_object = tracking->objects.first;
		     tracking_object;
		     tracking_object = tracking_object->next)
		{
			ListBase *tracksbase = BKE_tracking_object_get_tracks(tracking, tracking_object);
			MovieTrackingTrack *track;
			float reconstructed_camera_mat[4][4],
			      reconstructed_camera_imat[4][4];
			float (*vertex_obmat)[4];

			if ((tracking_object->flag & TRACKING_OBJECT_CAMERA) == 0) {
				BKE_tracking_camera_get_reconstructed_interpolate(tracking, tracking_object,
				                                                  CFRA, reconstructed_camera_mat);

				invert_m4_m4(reconstructed_camera_imat, reconstructed_camera_mat);
			}

			for (track = tracksbase->first; track; track = track->next) {
				float bundle_pos[3];

				if ((track->flag & TRACK_HAS_BUNDLE) == 0) {
					continue;
				}

				copy_v3_v3(bundle_pos, track->bundle_pos);
				if (tracking_object->flag & TRACKING_OBJECT_CAMERA) {
					vertex_obmat = orig_camera_mat;
				}
				else {
					mul_m4_v3(reconstructed_camera_imat, bundle_pos);
					vertex_obmat = obmat;
				}

				/* Use local values */
				mul_m4_v3(vertex_obmat, bundle_pos);
				retval |= snap_point_v3(
				        snpdt->mval, bundle_pos,
				        snpdt->pmat, snpdt->win_half, is_persp,
				        rv3d->clip, 4,
				        dist_2d, r_loc);
			}
		}

		if (retval) {
			*dist_px = sqrt(dist_2d[0] * dist_2d[0] + dist_2d[1] * dist_2d[1]);
			*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);
			return true;
		}
	}
	return retval;
}

static bool snapDerivedMesh(
        SnapObjectContext *sctx, SnapData *snpdt,
        SnapObjectBase *snap_obj, const unsigned int ob_index,
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float r_no[3], int *r_index,
        ListBase *r_hit_list)
{
	bool retval = false;

	DerivedMesh *dm = object_dm_final_get(sctx->scene, snap_obj->ob);

	if (dm->getNumVerts(dm) == 0) {
		return retval;
	}

	BoundBox *bb = BKE_object_boundbox_get(snap_obj->ob);

	if (snpdt->snap_to_flag & SCE_SELECT_FACE) {
		if (bb && !snp_snap_boundbox_raycast_test(snpdt, bb, snap_obj->mat)) {
			return retval;
		}
	}
	else {
		/* In vertex and edges you need to get the pixel distance from mval to BoundBox, see T46816 */
		if (bb && !snap_boundbox_nearest_test(snpdt, bb, snap_obj->mat, *dist_px)) {
			return retval;
		}
	}

	SnapDataMesh *snapdata;
	if (!snap_obj->data) {
		snapdata = snap_obj->data = BLI_memarena_calloc(sctx->cache.mem_arena, sizeof(*snapdata));
		snapdata->sd.type = SNAP_MESH;
		invert_m4_m4(snapdata->imat, snap_obj->mat);
		snapdata->has_looptris = true;
	}
	else {
		snapdata = snap_obj->data;
	}

	BVHTreeFromMesh *treedata = &snapdata->treedata;

	if ((snpdt->snap_to_flag & SCE_SELECT_FACE) && snapdata->has_looptris) {

		/* Adds data to cache */

		/* For any snap_to, the BVHTree of looptris will always be used */
		/* the tree is owned by the DM and may have been freed since we last used! */
		if (treedata->cached && !bvhcache_has_tree(dm->bvhCache, treedata->tree)) {
			free_bvhtree_from_mesh(treedata);
		}

		if (treedata->tree == NULL) {
			bvhtree_from_mesh_looptri(treedata, dm, 0.0f, 4, 6);

			/* Get mpoly array To check looptris reference (`DM_get_looptri_array`) */
			snapdata->mpoly = DM_get_poly_array(dm, &snapdata->poly_allocated);

			snapdata->has_looptris = treedata->tree != NULL;
		}
		else {
			/* The reference of the arrays can be lost in dm->release(dm); */
			if (!treedata->vert_allocated) {
				treedata->vert = DM_get_vert_array(dm, &treedata->vert_allocated);
			}
			if (!treedata->loop_allocated) {
				treedata->loop = DM_get_loop_array(dm, &treedata->loop_allocated);
			}
			if (!snapdata->poly_allocated) {
				snapdata->mpoly = DM_get_poly_array(dm, &snapdata->poly_allocated);
			}
			if (!treedata->looptri_allocated) {
				treedata->looptri = DM_get_looptri_array(
				        dm, treedata->vert,
				        snapdata->mpoly, dm->getNumPolys(dm),
				        treedata->loop, dm->getNumLoops(dm),
				        &treedata->looptri_allocated);
			}
		}

		float ray_start_local[3], ray_dir_local[3], local_scale, local_depth, len_diff = 0.0f;

		mul_v3_m4v3(ray_start_local, snap_obj->mat, snpdt->ray_start);
		mul_v3_mat3_m4v3(ray_dir_local, snapdata->imat, snpdt->ray_dir);

		/* local scale in normal direction */
		local_scale = normalize_v3(ray_dir_local);
		local_depth = *ray_depth;
		if (local_depth != BVH_RAYCAST_DIST_MAX) {
			local_depth *= local_scale;
		}

		/* Only use closer ray_start in case of ortho view! In perspective one, ray_start may already
		 * been *inside* boundbox, leading to snap failures (see T38409).
		 * Note also ar might be null (see T38435), in this case we assume ray_start is ok!
		 */
		if (snpdt->view_proj == VIEW_PROJ_ORTHO) {  /* do_ray_start_correction */
			if (bb) {
				len_diff = dist_aabb_to_plane(
				        bb->vec[0], bb->vec[6],
				        ray_start_local, ray_dir_local);

				if (len_diff < 0) len_diff = 0.0f;
			}
			else {
				/* TODO (Germano): Do with root */
			}

			/* You need to make sure that ray_start is really far away,
			 * because even in the Orthografic view, in some cases,
			 * the ray can start inside the object (see T50486) */
			if (len_diff > 400.0f) {
				/* We pass a temp ray_start, set from object's boundbox, to avoid precision issues with
				 * very far away ray_start values (as returned in case of ortho view3d), see T38358.
				 */
				float ray_org_local[3];
				copy_v3_v3(ray_org_local, snpdt->ray_origin);
				mul_m4_v3(snapdata->imat, ray_org_local);

				len_diff -= local_scale; /* make temp start point a bit away from bbox hit point. */
				madd_v3_v3v3fl(
				        ray_start_local, ray_org_local, ray_dir_local,
				        len_diff + snpdt->depth_range[0] * local_scale);
				local_depth -= len_diff;
			}
			else len_diff = 0.0f;
		}
		if (r_hit_list) {
			float timat[3][3];
			transpose_m3_m4(timat, snapdata->imat);

			struct RayCastAll_Data data;

			data.bvhdata = treedata;
			data.raycast_callback = treedata->raycast_callback;
			data.obmat = snap_obj->mat;
			data.timat = timat;
			data.len_diff = len_diff;
			data.local_scale = local_scale;
			data.ob = snap_obj->ob;
			data.ob_uuid = ob_index;
			data.hit_list = r_hit_list;
			data.retval = retval;

			BLI_bvhtree_ray_cast_all(
			        treedata->tree,
			        ray_start_local, ray_dir_local,
			        0.0f, *ray_depth, raycast_all_cb, &data);

			retval = data.retval;
		}
		else {
			BVHTreeRayHit hit = {.index = -1, .dist = local_depth};

			if (BLI_bvhtree_ray_cast(
			        treedata->tree,
			        ray_start_local, ray_dir_local,
			        0.0f, &hit, treedata->raycast_callback, treedata) != -1)
			{
				hit.dist += len_diff;
				hit.dist /= local_scale;
				if (hit.dist <= *ray_depth) {
					*ray_depth = hit.dist;
					copy_v3_v3(r_loc, hit.co);

					/* back to worldspace */
					mul_m4_v3(snap_obj->mat, r_loc);

					if (r_no) {
						float timat[3][3];
						transpose_m3_m4(timat, snapdata->imat);
						copy_v3_v3(r_no, hit.no);
						mul_m3_v3(timat, r_no);
						normalize_v3(r_no);
					}

					retval = true;

					if (r_index) {
						*r_index = treedata->looptri[hit.index].poly;
					}
				}
			}
		}
	}
	else { /* TODO (Germano): separate raycast and nearest */
		if (!treedata->vert_allocated) {
			treedata->vert = DM_get_vert_array(dm, &treedata->vert_allocated);
		}
		const MVert *mvert = treedata->vert;
		SnapPorjectVert *proj_vert = snapdata->sd.proj_vert;

		if (proj_vert == NULL) {
			bool is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
			float local_pmat[4][4];
			mul_m4_m4m4(local_pmat, snpdt->pmat, snap_obj->mat);
			/* The reference of the arrays can be lost in dm->release(dm); */
			int i, numVerts = dm->getNumVerts(dm);
			proj_vert = snapdata->sd.proj_vert = MEM_mallocN(sizeof(*proj_vert) * numVerts, __func__);
			const MVert *v = mvert;
			for (i = 0; i < numVerts; i++, v++) {
				snp_calc_project_vert(
				        NULL, 0, local_pmat, snpdt->win_half,
				        is_persp, v->co, &proj_vert[i]);
			}
		}

		short dist[2] = {*dist_px, *dist_px};
		int tmp_index = -1;

		if (snpdt->snap_to_flag & SCE_SELECT_VERTEX) {
			int i, numVerts = dm->getNumVerts(dm);
			SnapPorjectVert *v2d = proj_vert;
			for (i = 0; i < numVerts; i++, v2d++) {
				if (!v2d->is_visible) {
					continue;
				}

				short tmp_dist[2] = {
				    abs(snpdt->mval[0] - v2d->co[0]),
				    abs(snpdt->mval[1] - v2d->co[1]),
				};

				if (tmp_dist[0] <= dist[0] && tmp_dist[1] <= dist[1]) {
					copy_v2_v2_short(dist, tmp_dist);
					tmp_index = i;
				}
			}
			if (tmp_index != -1) {
				copy_v3_v3(r_loc, mvert[tmp_index].co);
				mul_m4_v3(snap_obj->mat, r_loc);

				if (r_no) {
					normal_short_to_float_v3(r_no, mvert[tmp_index].no);
					mul_transposed_mat3_m4_v3(snapdata->imat, r_no);
					normalize_v3(r_no);
				}

				*dist_px = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
				*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);

				retval = true;
			}
		}

		if (snpdt->snap_to_flag & SCE_SELECT_EDGE) {
			if (!treedata->edge_allocated) {
				treedata->edge = DM_get_edge_array(dm, &treedata->edge_allocated);
			}
			int i, numEdges = dm->getNumEdges(dm);
			float r_lambda;
			const MEdge *me = treedata->edge;
			for (i = 0; i < numEdges; i++, me++) {
				SnapPorjectVert *va2d = &proj_vert[me->v1];
				SnapPorjectVert *vb2d = &proj_vert[me->v2];

				if (!va2d->is_visible && !vb2d->is_visible) {
					continue;
				}

				if (snap_segment_v2v2(
				        snpdt->snap_to_flag, snpdt->mval, va2d, vb2d,
				        dist, &r_lambda))
				{
					tmp_index = i;
				}
			}
			if (tmp_index != -1) {
				me = &treedata->edge[tmp_index];
				float v1[3], v2[3];
				copy_v3_v3(v1, mvert[me->v1].co);
				copy_v3_v3(v2, mvert[me->v2].co);

				interp_v3_v3v3(r_loc, v1, v2, r_lambda);
				mul_m4_v3(snap_obj->mat, r_loc);

				if (r_no) {
					sub_v3_v3v3(r_no, v2, v1);
					mul_transposed_mat3_m4_v3(snapdata->imat, r_no);
					normalize_v3(r_no);
				}

				*dist_px = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
				*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);

				retval = true;
			}
		}
	}

	dm->release(dm);

	return retval;
}

static bool snapEditMesh(
        SnapObjectContext *sctx, SnapData *snpdt,
        SnapObjectBase *snap_obj, const unsigned int ob_index,
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float r_no[3], int *r_index,
        ListBase *r_hit_list)
{
	bool retval = false;

	SnapDataEditMesh *snapdata;
	if (!snap_obj->data) {
		snapdata = snap_obj->data = BLI_memarena_calloc(sctx->cache.mem_arena, sizeof(*snapdata));
		snapdata->sd.type = SNAP_MESH;
		snapdata->em = BKE_editmesh_from_object(snap_obj->ob);
		invert_m4_m4(snapdata->imat, snap_obj->mat);
	}
	else {
		snapdata = snap_obj->data;
	}

	if (snapdata->em->bm->totvert == 0) {
		return retval;
	}

	if (snpdt->snap_to_flag & SCE_SELECT_FACE) {

		if (snapdata->em->tottri == 0) {
			return retval;
		}

		BVHTreeFromEditMesh *treedata = &snapdata->treedata;
		if (treedata->tree == NULL) {
			int looptri_num_active = -1;
			BLI_bitmap *face_mask = NULL;
			if (sctx->callbacks.edit_mesh.test_face_fn) {
				face_mask = BLI_BITMAP_NEW(snapdata->em->tottri, __func__);
				looptri_num_active = BM_iter_mesh_bitmap_from_filter_tessface(
				        snapdata->em->bm, face_mask,
				        sctx->callbacks.edit_mesh.test_face_fn, sctx->callbacks.edit_mesh.user_data);
			}
			bvhtree_from_editmesh_looptri_ex(
			        treedata, snapdata->em, face_mask, looptri_num_active, 0.0f, 4, 6, NULL);

			if (face_mask) {
				MEM_freeN(face_mask);
			}
		}

		float ray_start_local[3], ray_dir_local[3], local_scale, local_depth, len_diff = 0.0f;

		mul_v3_m4v3(ray_start_local, snap_obj->mat, snpdt->ray_start);
		mul_v3_mat3_m4v3(ray_dir_local, snapdata->imat, snpdt->ray_dir);

		/* local scale in normal direction */
		local_scale = normalize_v3(ray_dir_local);
		local_depth = *ray_depth;
		if (local_depth != BVH_RAYCAST_DIST_MAX) {
			local_depth *= local_scale;
		}

		/* Only use closer ray_start in case of ortho view! In perspective one, ray_start
		 * may already been *inside* boundbox, leading to snap failures (see T38409).
		 * Note also ar might be null (see T38435), in this case we assume ray_start is ok!
		 */
		if (snpdt->view_proj == VIEW_PROJ_ORTHO) {  /* do_ray_start_correction */
			/* We *need* a reasonably valid len_diff in this case.
			 * Use BHVTree to find the closest face from ray_start_local.
			 */
			BVHTreeNearest nearest;
			nearest.index = -1;
			nearest.dist_sq = FLT_MAX;
			/* Compute and store result. */
			if (BLI_bvhtree_find_nearest(
			        treedata->tree, ray_start_local, &nearest, NULL, NULL) != -1)
			{
				float dvec[3];
				sub_v3_v3v3(dvec, nearest.co, ray_dir_local);
				len_diff = dot_v3v3(dvec, ray_dir_local);
				/* You need to make sure that ray_start is really far away,
				 * because even in the Orthografic view, in some cases,
				 * the ray can start inside the object (see T50486) */
				if (len_diff > 400.0f) {
					float ray_org_local[3];

					copy_v3_v3(ray_org_local, snpdt->ray_origin);
					mul_m4_v3(snapdata->imat, ray_org_local);

					/* We pass a temp ray_start, set from object's boundbox,
					 * to avoid precision issues with very far away ray_start values
					 * (as returned in case of ortho view3d), see T38358.
					 */
					len_diff -= local_scale; /* make temp start point a bit away from bbox hit point. */
					madd_v3_v3v3fl(
					        ray_start_local, ray_org_local, ray_dir_local,
					        len_diff + snpdt->depth_range[0] * local_scale);
					local_depth -= len_diff;
				}
				else len_diff = 0.0f;
			}
		}
		if (r_hit_list) {
			float timat[3][3];
			transpose_m3_m4(timat, snapdata->imat);

			struct RayCastAll_Data data;

			data.bvhdata = treedata;
			data.raycast_callback = treedata->raycast_callback;
			data.obmat = snap_obj->mat;
			data.timat = timat;
			data.len_diff = len_diff;
			data.local_scale = local_scale;
			data.ob = snap_obj->ob;
			data.ob_uuid = ob_index;
			data.hit_list = r_hit_list;
			data.retval = retval;

			BLI_bvhtree_ray_cast_all(
			        treedata->tree, ray_start_local, ray_dir_local, 0.0f,
			        *ray_depth, raycast_all_cb, &data);

			retval = data.retval;
		}
		else {
			BVHTreeRayHit hit = {.index = -1, .dist = local_depth};

			if (treedata->tree && BLI_bvhtree_ray_cast(
			        treedata->tree, ray_start_local, ray_dir_local, 0.0f,
			        &hit, treedata->raycast_callback, treedata) != -1)
			{
				hit.dist += len_diff;
				hit.dist /= local_scale;
				if (hit.dist <= *ray_depth) {
					*ray_depth = hit.dist;
					copy_v3_v3(r_loc, hit.co);

					/* back to worldspace */
					mul_m4_v3(snap_obj->mat, r_loc);

					if (r_no) {
						float timat[3][3];
						transpose_m3_m4(timat, snapdata->imat);

						copy_v3_v3(r_no, hit.no);
						mul_m3_v3(timat, r_no);
						normalize_v3(r_no);
					}

					retval = true;

					if (r_index) {
						*r_index = BM_elem_index_get(snapdata->em->looptris[hit.index][0]->f);
					}
				}
			}
		}
	}
	else { /* TODO (Germano): separate raycast from nearest */
		SnapPorjectVert *proj_vert = snapdata->sd.proj_vert;
		if (proj_vert == NULL) {
			bool is_persp = snpdt->view_proj == VIEW_PROJ_PERSP;
			float local_pmat[4][4];
			mul_m4_m4m4(local_pmat, snpdt->pmat, snap_obj->mat);
			/* The reference of the arrays can be lost in dm->release(dm); */
			int i, numVerts = snapdata->em->bm->totvert;
			proj_vert = snapdata->sd.proj_vert = MEM_mallocN(sizeof(*proj_vert) * numVerts, __func__);
			for (i = 0; i < numVerts; i++) {
				BMVert *v = BM_vert_at_index(snapdata->em->bm, i);
				if (BM_elem_flag_test(v, BM_ELEM_SELECT | BM_ELEM_HIDDEN)) {
					proj_vert[i].is_visible = false;
					continue;
				}
				snp_calc_project_vert(
				        NULL, 0, local_pmat, snpdt->win_half,
				        is_persp, v->co, &proj_vert[i]);
			}
		}

		short dist[2] = {*dist_px, *dist_px};
		int tmp_index = -1;

		if (snpdt->snap_to_flag & SCE_SELECT_VERTEX) {
			int i, numVerts = snapdata->em->bm->totvert;
			SnapPorjectVert *v2d = proj_vert;
			for (i = 0; i < numVerts; i++, v2d++) {
				if (!v2d->is_visible) {
					continue;
				}

				short tmp_dist[2] = {
					abs(snpdt->mval[0] - v2d->co[0]),
					abs(snpdt->mval[1] - v2d->co[1]),
				};

				if (tmp_dist[0] <= dist[0] && tmp_dist[1] <= dist[1]) {
					copy_v2_v2_short(dist, tmp_dist);
					tmp_index = i;
				}
			}
			if (tmp_index != -1) {
				BMVert *v = BM_vert_at_index(snapdata->em->bm, tmp_index);
				copy_v3_v3(r_loc, v->co);
				mul_m4_v3(snap_obj->mat, r_loc);

				if (r_no) {
					copy_v3_v3(r_no, v->no);
					mul_transposed_mat3_m4_v3(snapdata->imat, r_no);
					normalize_v3(r_no);
				}

				*dist_px = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
				*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);

				retval = true;
			}
		}

		if (snpdt->snap_to_flag & SCE_SELECT_EDGE) {
			int i, numEdges = snapdata->em->bm->totedge;
			float r_lambda;
			for (i = 0; i < numEdges; i++) {
				BMEdge *e = BM_edge_at_index(snapdata->em->bm, i);
				SnapPorjectVert *va2d = &proj_vert[BM_elem_index_get(e->v1)];
				SnapPorjectVert *vb2d = &proj_vert[BM_elem_index_get(e->v2)];

				if (!va2d->is_visible && !vb2d->is_visible) {
					continue;
				}

				if (snap_segment_v2v2(
				        snpdt->snap_to_flag, snpdt->mval, va2d, vb2d,
				        dist, &r_lambda))
				{
					tmp_index = i;
				}
			}

			if (tmp_index != -1) {
				BMEdge *e = BM_edge_at_index(snapdata->em->bm, tmp_index);
				float v1[3], v2[3];
				copy_v3_v3(v1, e->v1->co);
				copy_v3_v3(v2, e->v2->co);

				interp_v3_v3v3(r_loc, v1, v2, r_lambda);
				mul_m4_v3(snap_obj->mat, r_loc);

				if (r_no) {
					sub_v3_v3v3(r_no, v2, v1);
					mul_transposed_mat3_m4_v3(snapdata->imat, r_no);
					normalize_v3(r_no);
				}
				*dist_px = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
				*ray_depth = depth_get(r_loc, snpdt->ray_start, snpdt->ray_dir);

				retval = true;
			}
		}
	}

	return retval;
}

/**
 * \param use_obedit: Uses the coordinates of BMesh (if any) to do the snapping;
 *
 * \note Duplicate args here are documented at #snapObjectsRay
 */
static bool snapObject(
        SnapObjectContext *sctx, SnapData *snpdt,
        SnapObjectBase *snap_obj, const unsigned int ob_index,
        bool use_obedit,
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float r_no[3], int *r_index,
        SnapObjectBase **r_ob, ListBase *r_hit_list)
{
	bool retval = false;

	if (snpdt->test_occlusion &&
		snpdt->snap_to_flag == SCE_SELECT_FACE &&
		(sctx->v3d_data.v3d->drawtype < OB_SOLID ||
		snap_obj->ob->dt < OB_SOLID))
	{
		return retval;
	}

	switch (snap_obj->ob->type) {
		case OB_MESH:
		{
			if (use_obedit) {
				retval = snapEditMesh(
				        sctx, snpdt, snap_obj, ob_index,
				        ray_depth, dist_px,
				        r_loc, r_no, r_index,
				        r_hit_list);
			}
			else {
				retval = snapDerivedMesh(
				        sctx, snpdt, snap_obj, ob_index,
				        ray_depth, dist_px,
				        r_loc, r_no,
				        r_index, r_hit_list);
			}
			break;
		}
		case OB_ARMATURE:
			retval = snapArmature(
			        sctx, snpdt,
			        snap_obj->ob, snap_obj->ob->data, snap_obj->mat,
			        ray_depth, dist_px,
			        r_loc, r_no);
			break;
		case OB_CURVE:
			retval = snapCurve(
			        sctx, snpdt,
			        snap_obj->ob, snap_obj->ob->data, snap_obj->mat,
			        ray_depth, dist_px,
			        r_loc, r_no);
			break;
		case OB_EMPTY:
			retval = snapEmpty(
			        sctx, snpdt,
			        snap_obj->ob, snap_obj->mat,
			        ray_depth, dist_px,
			        r_loc, r_no);
			break;
		case OB_CAMERA:
			retval = snapCamera(
			        sctx, snpdt, snap_obj->ob, snap_obj->mat,
			        ray_depth, dist_px,
			        r_loc, r_no);
			break;
	}

	if (retval) {
		if (r_ob) {
			*r_ob = snap_obj;
		}
	}

	return retval;
}

/**
 * Main Snapping Function
 * ======================
 *
 * Walks through all objects in the scene to find the closest snap element ray.
 *
 * \param sctx: Snap context to store data.
 * \param snpdt: struct generated in `snapdata_init_ray/v3d`.
 * \param snap_select: from enum SnapSelect.
 * \param use_object_edit_cage: Uses the coordinates of BMesh (if any) to do the snapping.
 *
 * Read/Write Args
 * ---------------
 *
 * \param ray_depth: maximum depth allowed for r_co, elements deeper than this value will be ignored.
 * \param dist_px: Maximum threshold distance (in pixels).
 *
 * Output Args
 * -----------
 *
 * \param r_loc: Hit location.
 * \param r_no: Hit normal (optional).
 * \param r_index: Hit index or -1 when no valid index is found.
 * (currently only set to the polygon index when when using ``snap_to == SCE_SNAP_MODE_FACE``).
 * \param r_ob: Hit object.
 * \param r_obmat: Object matrix (may not be #Object.obmat with dupli-instances).
 * \param r_hit_list: List of #SnapObjectHitDepth (caller must free).
 *
 */
static bool snapObjectsRay(
        SnapObjectContext *sctx, SnapData *snpdt,
        /* read/write args */
        float *ray_depth, float *dist_px,
        /* return args */
        float r_loc[3], float r_no[3], int *r_index,
        SnapObjectBase **r_ob, ListBase *r_hit_list)
{
	bool retval = false;

	unsigned int ob_index = 0;
	if (!BLI_listbase_is_empty(&sctx->cache.object_map)) {
		for (SnapObjectBase *snap_obj = sctx->cache.object_map.first; snap_obj; snap_obj = snap_obj->next) {
			retval |= snapObject(
			        sctx, snpdt, snap_obj, ob_index++, snap_obj->is_obedit,
			        ray_depth, dist_px,
			        r_loc, r_no, r_index, r_ob, r_hit_list);
		}
	}

	return retval;
}

/**
 *
 */
static void create_object_scene_list(
        ListBase *obj_listmap,
        SnapObjectContext *sctx,
        const SnapSelect snap_select,
        const bool use_object_edit_cage)
{
	bool retval = false;

	unsigned int ob_index = 0;
	Object *obedit = use_object_edit_cage ? sctx->scene->obedit : NULL;

	/* Need an exception for particle edit because the base is flagged with BA_HAS_RECALC_DATA
	 * which makes the loop skip it, even the derived mesh will never change
	 *
	 * To solve that problem, we do it first as an exception.
	 * */
	Base *base_act = sctx->scene->basact;
	if (base_act && base_act->object && base_act->object->mode & OB_MODE_PARTICLE_EDIT) {
		Object *ob = base_act->object;

		SnapObjectBase *snap_obj = MEM_mallocN(sizeof(*snap_obj), __func__);

		snap_obj->ob = ob;
		snap_obj->is_obedit = false;
		copy_m4_m4(snap_obj->mat, ob->obmat);
		snap_obj->data = NULL;

		BLI_addtail(obj_listmap, snap_obj);
	}

	bool ignore_object_selected = false, ignore_object_active = false;
	switch (snap_select) {
		case SNAP_ALL:
			break;
		case SNAP_NOT_SELECTED:
			ignore_object_selected = true;
			break;
		case SNAP_NOT_ACTIVE:
			ignore_object_active = true;
			break;
	}
	for (Base *base = sctx->scene->base.first; base != NULL; base = base->next) {
		if ((BASE_VISIBLE_BGMODE(sctx->v3d_data.v3d, sctx->scene, base)) &&
		    (base->flag & (BA_HAS_RECALC_OB | BA_HAS_RECALC_DATA)) == 0 &&

		    !((ignore_object_selected && (base->flag & (SELECT | BA_WAS_SEL))) ||
		      (ignore_object_active && base == base_act)))
		{
			Object *ob = base->object;

			if (ob->transflag & OB_DUPLI) {
				DupliObject *dupli_ob;
				ListBase *lb = object_duplilist(sctx->bmain->eval_ctx, sctx->scene, ob);

				for (dupli_ob = lb->first; dupli_ob; dupli_ob = dupli_ob->next) {
					bool use_obedit_dupli = (obedit && dupli_ob->ob->data == obedit->data);
					Object *obj_dupli = (use_obedit_dupli) ? obedit : dupli_ob->ob;

					SnapObjectBase *snap_obj = MEM_mallocN(sizeof(*snap_obj), __func__);

					snap_obj->ob = obj_dupli;
					snap_obj->is_obedit = use_obedit_dupli;
					copy_m4_m4(snap_obj->mat, dupli_ob->mat);
					snap_obj->data = NULL;

					BLI_addtail(obj_listmap, snap_obj);
				}

				free_object_duplilist(lb);
			}

			bool use_obedit = (obedit != NULL) && (ob->data == obedit->data);
			Object *ob_snap = use_obedit ? obedit : ob;

			SnapObjectBase *snap_obj = MEM_mallocN(sizeof(*snap_obj), __func__);

			snap_obj->ob = ob_snap;
			snap_obj->is_obedit = use_obedit;
			copy_m4_m4(snap_obj->mat, ob->obmat);
			snap_obj->data = NULL;

			BLI_addtail(obj_listmap, snap_obj);
		}
	}
}

/** \} */


/* -------------------------------------------------------------------- */

/** \name Public Object Snapping API
 * \{ */

SnapObjectContext *ED_transform_snap_object_context_create(
        const struct SnapObjectParams *params,
        Main *bmain, Scene *scene, int flag)
{
	SnapObjectContext *sctx = MEM_callocN(sizeof(*sctx), __func__);

	sctx->flag = flag;

	sctx->bmain = bmain;
	sctx->scene = scene;

	sctx->cache.mem_arena = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, __func__);
	create_object_scene_list(&sctx->cache.object_map, sctx, params->snap_select, params->use_object_edit_cage);

	return sctx;
}

SnapObjectContext *ED_transform_snap_object_context_create_view3d(
        const struct SnapObjectParams *params,
        Main *bmain, Scene *scene, int flag,
        /* extra args for view3d */
        const ARegion *ar, const View3D *v3d)
{
	SnapObjectContext *sctx = ED_transform_snap_object_context_create(params, bmain, scene, flag);

	sctx->use_v3d = true;
	sctx->v3d_data.ar = ar;
	sctx->v3d_data.v3d = v3d;

	return sctx;
}

static void snap_object_data_free(void *sod_v)
{
	SnapObjectData *sod = sod_v;
	if (sod->proj_vert) {
		MEM_freeN(sod->proj_vert);
	}
	switch (sod->type) {
		case SNAP_MESH:
		{
			SnapDataMesh *sodtype = sod_v;
			if (sodtype->treedata.tree) {
				free_bvhtree_from_mesh(&sodtype->treedata);
				if (sodtype->poly_allocated) {
					MEM_freeN(sodtype->mpoly);
				}
			}
			break;
		}
		case SNAP_EDIT_MESH:
		{
			SnapDataEditMesh *sodtype = sod_v;
			if (sodtype->treedata.tree) {
				free_bvhtree_from_editmesh(&sodtype->treedata);
			}
			break;
		}
	}
}

void ED_transform_snap_object_context_destroy(SnapObjectContext *sctx)
{
	BLI_freelistN(&sctx->cache.object_map);
	BLI_memarena_free(sctx->cache.mem_arena);

	MEM_freeN(sctx);
}

void ED_transform_snap_object_context_set_editmesh_callbacks(
        SnapObjectContext *sctx,
        bool (*test_vert_fn)(BMVert *, void *user_data),
        bool (*test_edge_fn)(BMEdge *, void *user_data),
        bool (*test_face_fn)(BMFace *, void *user_data),
        void *user_data)
{
	sctx->callbacks.edit_mesh.test_vert_fn = test_vert_fn;
	sctx->callbacks.edit_mesh.test_edge_fn = test_edge_fn;
	sctx->callbacks.edit_mesh.test_face_fn = test_face_fn;

	sctx->callbacks.edit_mesh.user_data = user_data;
}

bool ED_transform_snap_object_project_ray_ex(
        SnapObjectContext *sctx,
        const float ray_start[3], const float ray_normal[3],
        float *ray_depth,
        float r_loc[3], float r_no[3], int *r_index,
        Object **r_ob, float r_obmat[4][4])
{
	SnapData snpdt;
	snapdata_init_ray(&snpdt, ray_start, ray_normal);

	SnapObjectBase *snap_obj = NULL;

	if (snapObjectsRay(
	        sctx, &snpdt,
	        ray_depth, NULL,
	        r_loc, r_no, r_index, &snap_obj, NULL))
	{
		if (r_index) {
			/* Restore index exposed by polys in in bpy.object.data */
			SnapObjectData *sod = snap_obj->data;
			if (sod->type == SNAP_MESH) {

				DerivedMesh *dm = object_dm_final_get(sctx->scene, snap_obj->ob);

				const int *index_mp_to_orig = dm->getPolyDataArray(dm, CD_ORIGINDEX);
				*r_index = index_mp_to_orig ? index_mp_to_orig[*r_index] : *r_index;

				dm->release(dm);
			}
		}
		if (r_ob) {
			*r_ob = snap_obj->ob;
		}

		return true;
	}

	return false;
}

/**
 * Fill in a list of all hits.
 *
 * \param ray_depth: Only depths in this range are considered, -1.0 for maximum.
 * \param sort: Optionally sort the hits by depth.
 * \param r_hit_list: List of #SnapObjectHitDepth (caller must free).
 */
bool ED_transform_snap_object_project_ray_all(
        SnapObjectContext *sctx,
        const float ray_start[3], const float ray_normal[3],
        float ray_depth, bool sort,
        ListBase *r_hit_list)
{
	const float depth_range[2] = {0.0f, FLT_MAX};
	if (ray_depth == -1.0f) {
		ray_depth = BVH_RAYCAST_DIST_MAX;
	}

#ifdef DEBUG
	float ray_depth_prev = ray_depth;
#endif

	SnapData snpdt;
	snapdata_init_ray(&snpdt, ray_start, ray_normal);

	bool retval = snapObjectsRay(
	        sctx, &snpdt, &ray_depth, NULL,
	        NULL, NULL, NULL, NULL, r_hit_list);

	/* meant to be readonly for 'all' hits, ensure it is */
#ifdef DEBUG
	BLI_assert(ray_depth_prev == ray_depth);
#endif

	if (sort) {
		BLI_listbase_sort(r_hit_list, hit_depth_cmp_cb);
	}

	return retval;
}

/**
 * Convenience function for snap ray-casting.
 *
 * Given a ray, cast it into the scene (snapping to faces).
 *
 * \return Snap success
 */
static bool transform_snap_context_project_ray_impl(
        SnapObjectContext *sctx,
        const float ray_start[3], const float ray_normal[3], float *ray_depth,
        float r_co[3], float r_no[3])
{
	/* try snap edge, then face if it fails */
	bool ret = ED_transform_snap_object_project_ray_ex(
	        sctx,
	        ray_start, ray_normal, ray_depth,
	        r_co, r_no, NULL,
	        NULL, NULL);

	return ret;
}

bool ED_transform_snap_object_project_ray(
        SnapObjectContext *sctx,
        const float ray_origin[3], const float ray_direction[3], float *ray_depth,
        float r_co[3], float r_no[3])
{
	float ray_depth_fallback;
	if (ray_depth == NULL) {
		ray_depth_fallback = BVH_RAYCAST_DIST_MAX;
		ray_depth = &ray_depth_fallback;
	}

	return transform_snap_context_project_ray_impl(
	        sctx,
	        ray_origin, ray_direction, ray_depth,
	        r_co, r_no);
}

static bool transform_snap_context_project_view3d_mixed_impl(
        SnapObjectContext *sctx,
        SnapData *snpdt, float *dist_px, float *ray_depth,
        float r_co[3], float r_no[3], int *r_index)
{
	BLI_assert(snpdt->snap_to_flag != 0);
	BLI_assert((snpdt->snap_to_flag & ~(1 | 2 | 4)) == 0);

	bool is_hit = false;
	int t_index;
	SnapObjectBase *snap_obj = NULL;
	float t_no[3];

	if (snpdt->test_occlusion || snpdt->snap_to_flag & SCE_SELECT_FACE) {
		unsigned short tmp_snap_to_flag = snpdt->snap_to_flag;
		snpdt->snap_to_flag = SCE_SELECT_FACE;
		if (snapObjectsRay(
		        sctx, snpdt,
		        ray_depth, dist_px,
		        r_co, t_no, &t_index, &snap_obj, NULL))
		{
			is_hit = (tmp_snap_to_flag & SCE_SELECT_FACE) != 0;

			/* Get new clip plane to simule occlusion */
			if (tmp_snap_to_flag & (SCE_SELECT_EDGE | SCE_SELECT_VERTEX)) {

				snpdt->dinamic_plane = MEM_mallocN(sizeof(*snpdt->dinamic_plane), __func__);
				float normal_local[3], plane_no[3], far_vert[3];

				copy_v3_v3(plane_no, t_no);
				if (dot_v3v3(plane_no, snpdt->ray_dir) > 0) {
					negate_v3(plane_no);
				}

				SnapObjectData *sod = snap_obj->data;
				if (sod->type == SNAP_MESH) {
					SnapDataMesh *snapdata = snap_obj->data;
					copy_v3_v3(normal_local, plane_no);
					mul_m4_v3(snapdata->imat, normal_local);

					DerivedMesh *dm = object_dm_final_get(sctx->scene, snap_obj->ob);

					/* The reference of the arrays can be lost in dm->release(dm); */
					if (!snapdata->treedata.vert_allocated) {
						snapdata->treedata.vert = DM_get_vert_array(dm, &snapdata->treedata.vert_allocated);
					}
					if (!snapdata->treedata.loop_allocated) {
						snapdata->treedata.loop = DM_get_loop_array(dm, &snapdata->treedata.loop_allocated);
					}
					if (!snapdata->poly_allocated) {
						snapdata->mpoly = DM_get_poly_array(dm, &snapdata->poly_allocated);
					}

					MPoly *mp = &snapdata->mpoly[t_index];

					int loopstart = mp->loopstart;
					int totloop = mp->totloop;

					copy_v3_v3(far_vert,
					           snapdata->treedata.vert[snapdata->treedata.loop[loopstart + totloop - 1].v].co);

					float iter_dist, far_dist = dot_v3v3(far_vert, normal_local);

					const MLoop *ml = &snapdata->treedata.loop[loopstart];
					for (int i = 1; i < totloop; i++, ml++) {
						iter_dist = dot_v3v3(snapdata->treedata.vert[ml->v].co, normal_local);
						if (iter_dist < far_dist) {
							far_dist = iter_dist;
							copy_v3_v3(far_vert, snapdata->treedata.vert[ml->v].co);
						}
					}

					dm->release(dm);
				}
				else { //if (sod->type == SNAP_EDIT_MESH)
					SnapDataEditMesh *snapdata = snap_obj->data;
					copy_v3_v3(normal_local, plane_no);
					mul_m4_v3(snapdata->imat, normal_local);

					BMEditMesh *em = snapdata->em;
					BMFace *f = BM_face_at_index(em->bm, t_index);

					BMLoop *l_iter, *l_first;
					l_first = BM_FACE_FIRST_LOOP(f);

					copy_v3_v3(far_vert, l_first->v->co);
					float iter_dist, far_dist = dot_v3v3(far_vert, normal_local);
					l_iter = l_first->next;

					do {
						iter_dist = dot_v3v3(l_iter->v->co, normal_local);
						if (iter_dist < far_dist) {
							far_dist = iter_dist;
							copy_v3_v3(far_vert, l_iter->v->co);
						}
					} while ((l_iter = l_iter->next) != l_first);
				}

				mul_m4_v3(snap_obj->mat, far_vert);

				plane_from_point_normal_v3(
				        snpdt->dinamic_plane[0], far_vert, plane_no);

				/* Slightly move the clip plane away since there was no snap in the polygon (TODO) */
				snpdt->dinamic_plane[0][3] += 0.000005;
			}
		}

		snpdt->snap_to_flag = tmp_snap_to_flag & ~SCE_SELECT_FACE;
	}


	if (snpdt->snap_to_flag) {
		BLI_assert(dist_px != NULL);
		if (snapObjectsRay(
			        sctx, snpdt,
			        ray_depth, dist_px,
			        r_co, t_no, &t_index, NULL, NULL))
		{
			is_hit = true;
		}
	}

	if (r_no) {
		copy_v3_v3(r_no, t_no);
	}
	if (r_index) {
		*r_index = t_index;
	}

	return is_hit;
}

/**
 * Convenience function for performing snapping.
 *
 * Given a 2D region value, snap to vert/edge/face.
 *
 * \param sctx: Snap context.
 * \param mval_fl: Screenspace coordinate.
 * \param dist_px: Maximum distance to snap (in pixels).
 * \param use_depth: Snap to the closest element, use when using more than one snap type.
 * \param r_co: hit location.
 * \param r_no: hit normal (optional).
 * \return Snap success
 */
bool ED_transform_snap_object_project_view3d_mixed(
        SnapObjectContext *sctx,
        const unsigned short snap_to_flag,
        const float mval_fl[2], float *dist_px,
        bool use_depth,
        float r_co[3], float r_no[3])
{
	float ray_depth = BVH_RAYCAST_DIST_MAX;

	SnapData snpdt;
	if (!snapdata_init_v3d(&snpdt, sctx, snap_to_flag, mval_fl, &ray_depth)) {
		return false;
	}

	snpdt.snap_to_flag = snap_to_flag;
	snpdt.test_occlusion = use_depth;

	bool ret = transform_snap_context_project_view3d_mixed_impl(
	        sctx, &snpdt, dist_px, &ray_depth,
	        r_co, r_no, NULL);

	if (snpdt.rv3d_clip) {
		MEM_freeN(snpdt.rv3d_clip);
	}
	if (snpdt.dinamic_plane) {
		MEM_freeN(snpdt.dinamic_plane);
	}

	return ret;
}

bool ED_transform_snap_object_project_view3d_ex(
        SnapObjectContext *sctx,
        const unsigned short snap_to,
        const float mval[2], float *dist_px,
        float *ray_depth,
        float r_loc[3], float r_no[3], int *r_index)
{
	unsigned short snap_to_flag;
	switch (snap_to) {
		case SCE_SNAP_MODE_FACE:
			snap_to_flag = SCE_SELECT_FACE;
			break;
		case SCE_SNAP_MODE_VERTEX:
			snap_to_flag = SCE_SELECT_VERTEX;
			break;
		case SCE_SNAP_MODE_EDGE:
			snap_to_flag = SCE_SELECT_EDGE;
			break;
		default:
			return false;
	}

	float ray_depth_fallback = BVH_RAYCAST_DIST_MAX;
	if (ray_depth == NULL) {
		ray_depth = &ray_depth_fallback;
	}

	SnapData snpdt;
	if (!snapdata_init_v3d(&snpdt, sctx, snap_to_flag, mval, ray_depth)) {
		return false;
	}

	if (snap_to_flag == SCE_SELECT_FACE) {
		/* Users may want to snap-to-face of objects that are in wireframe or BoundBox mode
		 * (this follows the previous behavior of blender) */
		snpdt.test_occlusion = false;
	}

	bool ret = transform_snap_context_project_view3d_mixed_impl(
	        sctx, &snpdt, dist_px, ray_depth,
	        r_loc, r_no, r_index);

	if (snpdt.rv3d_clip) {
		MEM_freeN(snpdt.rv3d_clip);
	}
	if (snpdt.dinamic_plane) {
		MEM_freeN(snpdt.dinamic_plane);
	}

	return ret;
}

bool ED_transform_snap_object_project_view3d(
        SnapObjectContext *sctx,
        const unsigned short snap_to,
        const float mval[2], float *dist_px,
        float *ray_depth,
        float r_loc[3], float r_no[3])
{
	return ED_transform_snap_object_project_view3d_ex(
	        sctx,
	        snap_to,
	        mval, dist_px,
	        ray_depth,
	        r_loc, r_no, NULL);
}

/**
 * see: #ED_transform_snap_object_project_ray_all
 */
bool ED_transform_snap_object_project_all_view3d_ex(
        SnapObjectContext *sctx,
        const float mval[2],
        float ray_depth, bool sort,
        ListBase *r_hit_list)
{
	float ray_start[3], ray_normal[3];

	if (!ED_view3d_win_to_ray_ex(
	        sctx->v3d_data.ar, sctx->v3d_data.v3d,
	        mval, NULL, ray_normal, ray_start, true))
	{
		return false;
	}

	return ED_transform_snap_object_project_ray_all(
	        sctx,
	        ray_start, ray_normal, ray_depth, sort,
	        r_hit_list);
}

/** \} */
