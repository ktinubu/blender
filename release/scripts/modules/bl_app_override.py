
# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8-80 compliant>

"""
Module to manage overriding various parts of Blender.

Intended for use with 'app_templates', though it can be used from anywhere.
"""

class class_match:
    # TODO, how to check these aren't from add-ons,
    # templates might need to un-register while filtering.

    @staticmethod
    def _rna_subclasses(cls_parent):
        for cls in cls_parent.__subclasses__():
            yield cls

    @staticmethod
    def panel(*, bl_category=None, bl_space_type=None, bl_region_type=None):
        # None or set
        import bpy
        for cls in class_match._rna_subclasses(bpy.types.Panel):
            if bl_category is not None and cls.bl_category not in bl_category:
                continue
            if bl_space_type is not None and cls.bl_space_type not in bl_space_type:
                continue
            if bl_region_type is not None and cls.bl_region_type not in bl_region_type:
                continue
            yield cls
