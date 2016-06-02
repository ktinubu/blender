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

/** \file blender/depsgraph/intern/build/deg_builder.cc
 *  \ingroup depsgraph
 */

#include "intern/builder/deg_builder.h"

// TODO(sergey): Use own wrapper over STD.
#include <stack>

#include "DNA_anim_types.h"

#include "BLI_utildefines.h"
#include "BLI_ghash.h"

#include "intern/depsgraph.h"
#include "intern/depsgraph_types.h"
#include "intern/nodes/deg_node.h"
#include "intern/nodes/deg_node_component.h"
#include "intern/nodes/deg_node_operation.h"

#include "util/deg_util_foreach.h"

namespace DEG {

string deg_fcurve_id_name(const FCurve *fcu)
{
	char index_buf[32];
	// TODO(sergey): Use int-to-string utility or so.
	BLI_snprintf(index_buf, sizeof(index_buf), "[%d]", fcu->array_index);
	return string(fcu->rna_path) + index_buf;
}

void deg_graph_build_finalize(Depsgraph *graph)
{
	std::stack<OperationDepsNode *> stack;

	foreach (OperationDepsNode *node, graph->operations) {
		node->done = 0;
		node->num_links_pending = 0;
		foreach (DepsRelation *rel, node->inlinks) {
			if ((rel->from->type == DEPSNODE_TYPE_OPERATION) &&
			    (rel->flag & DEPSREL_FLAG_CYCLIC) == 0)
			{
				++node->num_links_pending;
			}
		}
		if (node->num_links_pending == 0) {
			stack.push(node);
		}
		IDDepsNode *id_node = node->owner->owner;
		id_node->id->tag |= LIB_TAG_DOIT;
	}

	while (!stack.empty()) {
		OperationDepsNode *node = stack.top();
		if (node->done == 0 && node->outlinks.size() != 0) {
			foreach (DepsRelation *rel, node->outlinks) {
				if (rel->to->type == DEPSNODE_TYPE_OPERATION) {
					OperationDepsNode *to = (OperationDepsNode *)rel->to;
					if ((rel->flag & DEPSREL_FLAG_CYCLIC) == 0) {
						BLI_assert(to->num_links_pending > 0);
						--to->num_links_pending;
					}
					if (to->num_links_pending == 0) {
						stack.push(to);
					}
				}
			}
			node->done = 1;
		}
		else {
			stack.pop();
			IDDepsNode *id_node = node->owner->owner;
			foreach (DepsRelation *rel, node->outlinks) {
				if (rel->to->type == DEPSNODE_TYPE_OPERATION) {
					OperationDepsNode *to = (OperationDepsNode *)rel->to;
					IDDepsNode *id_to = to->owner->owner;
					id_node->layers |= id_to->layers;
				}
			}
		}
	}

	/* Re-tag IDs for update if it was tagged before the relations update tag. */
	GHASH_FOREACH_BEGIN(IDDepsNode *, id_node, graph->id_hash)
	{
		ID *id = id_node->id;
		if (id->tag & LIB_TAG_ID_RECALC_ALL &&
		    id->tag & LIB_TAG_DOIT)
		{
			id_node->tag_update(graph);
			id->tag &= ~LIB_TAG_DOIT;
		}
		id_node->finalize_build();
	}
	GHASH_FOREACH_END();
}

DEG::eDepsNode_Type deg_build_get_component_type(
        eDepsComponent component)
{
	switch (component) {
		case DEG_COMPONENT_PARAMETERS: return DEG::DEPSNODE_TYPE_PARAMETERS;
		case DEG_COMPONENT_PROXY: return DEG::DEPSNODE_TYPE_PROXY;
		case DEG_COMPONENT_ANIMATION: return DEG::DEPSNODE_TYPE_ANIMATION;
		case DEG_COMPONENT_TRANSFORM: return DEG::DEPSNODE_TYPE_TRANSFORM;
		case DEG_COMPONENT_GEOMETRY: return DEG::DEPSNODE_TYPE_GEOMETRY;
		case DEG_COMPONENT_SEQUENCER: return DEG::DEPSNODE_TYPE_SEQUENCER;
		case DEG_COMPONENT_EVAL_POSE: return DEG::DEPSNODE_TYPE_EVAL_POSE;
		case DEG_COMPONENT_BONE: return DEG::DEPSNODE_TYPE_BONE;
		case DEG_COMPONENT_EVAL_PARTICLES: return DEG::DEPSNODE_TYPE_EVAL_PARTICLES;
		case DEG_COMPONENT_SHADING: return DEG::DEPSNODE_TYPE_SHADING;
	}
	return DEG::DEPSNODE_TYPE_UNDEFINED;
}

}  // namespace DEG