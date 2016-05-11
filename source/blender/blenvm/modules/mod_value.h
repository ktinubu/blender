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
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __MOD_VALUE_H__
#define __MOD_VALUE_H__

#include "mod_defines.h"

BVM_MOD_NAMESPACE_BEGIN

BVM_MOD_FUNCTION("VALUE_FLOAT")
void VALUE_FLOAT(float &result, float value)
{
	result = value;
}

BVM_MOD_FUNCTION("VALUE_FLOAT3")
void VALUE_FLOAT3(float3 &result, const float3 &value)
{
	result = value;
}

BVM_MOD_FUNCTION("VALUE_FLOAT4")
void VALUE_FLOAT4(float4 &result, const float4 &value)
{
	result = value;
}

BVM_MOD_FUNCTION("VALUE_INT")
void VALUE_INT(int &result, int value)
{
	result = value;
}

BVM_MOD_FUNCTION("VALUE_MATRIX44")
void VALUE_MATRIX44(matrix44 &result, const matrix44 &value)
{
	result = value;
}

BVM_MOD_NAMESPACE_END

#endif /* __MOD_VALUE_H__ */