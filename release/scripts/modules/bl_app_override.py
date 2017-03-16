
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


# TODO, how to check these aren't from add-ons.
# templates might need to un-register while filtering.
def class_filter(cls_parent, **kw):
    kw = kw.copy()
    white_list = kw.pop("white_list", None)
    black_list = kw.pop("black_list", None)
    kw_items = tuple(kw.items())
    for cls in cls_parent.__subclasses__():
        if black_list is not None and cls.__name__ in black_list:
            continue
        if ((white_list is not None and cls.__name__ is white_list) or
                all((getattr(cls, attr) in expect) for attr, expect in kw_items)):
            yield cls
