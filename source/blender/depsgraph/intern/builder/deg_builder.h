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
 * The Original Code is Copyright (C) 2016 Blender Foundation.
 * All rights reserved.
 *
 * Original Author: Sergey Sharybin
 * Contributor(s): None Yet
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/depsgraph/intern/build/deg_builder.h
 *  \ingroup depsgraph
 */

#pragma once

extern "C" {
#include "DEG_depsgraph_build.h"
}

#include "intern/depsgraph_types.h"

struct FCurve;

namespace DEG {

struct Depsgraph;

/* Get unique identifier for FCurves and Drivers */
string deg_fcurve_id_name(const FCurve *fcu);

void deg_graph_build_finalize(struct Depsgraph *graph);

DEG::eDepsNode_Type deg_build_get_component_type(eDepsComponent component);

}  // namespace DEG