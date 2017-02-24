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

/** \file DNA_workspace_types.h
 *  \ingroup DNA
 *
 * Only use with API in BKE_workspace.h!
 */

#ifndef __DNA_WORKSPACE_TYPES_H__
#define __DNA_WORKSPACE_TYPES_H__

#if !defined(NAMESPACE_WORKSPACE) && !defined(NAMESPACE_DNA)
#  error "This file shouldn't be included outside of workspace namespace."
#endif

/**
 * \brief Wrapper for bScreen.
 *
 * bScreens are IDs and thus stored in a main list-base. We also want to store a list-base of them within the
 * workspace (so each workspace can have its own set of screen-layouts) which would mess with the next/prev pointers.
 * So we use this struct to wrap a bScreen pointer with another pair of next/prev pointers.
 *
 * We could also use LinkNode for this but in future we may want to move stuff from bScreen to this level.
 */
typedef struct WorkSpaceLayout {
	struct WorkSpaceLayout *next, *prev;

	struct WorkSpaceLayoutType *type;
	struct bScreen *screen;
} WorkSpaceLayout;

typedef struct WorkSpaceLayoutType {
	struct WorkSpaceLayoutType *next, *prev;
	const char *name;
} WorkSpaceLayoutType;

typedef struct WorkSpace {
	ID id;

	ListBase layout_types;
	ListBase layouts;
	struct WorkSpaceLayoutType *act_layout_type;
	struct WorkSpaceLayout *act_layout;
	/* We can't switch layout from within handlers since handler loop
	 * heavily depends on layout, so we store it here and change later */
	struct WorkSpaceLayout *new_layout;

	int object_mode; /* enum ObjectMode */
	int pad;
} WorkSpace;

/**
 * This struct is the bridge between workspaces and the entity type they belong to, currently wmWindow.
 * It makes it possible to manage workspace data completely on workspace level, totally separate from wmWindow.
 */
typedef struct WorkSpaceHook {
	WorkSpace *act_workspace;
	/* We can't switch workspace from within handlers since handler loop
	 * heavily depends on workspace, so we store it here and change later */
	WorkSpace *new_workspace; /* temporary when switching */

	/* To support opening a workspace in multiple windows while keeping the individual layouts independent, each
	 * window stores a list of layouts that is synced with a list of layout-type definitions from the workspace */
	ListBase layouts;
} WorkSpaceHook;

#endif /* __DNA_WORKSPACE_TYPES_H__ */
