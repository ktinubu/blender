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

class_store = []

import bl_app_override

def register():
    import bpy
    print("Template Register", __file__)

    class_store.clear()

    class_store.extend(bl_app_override.class_match.panel())
    unregister = bpy.utils.unregister_class
    for cls in class_store:
        unregister(cls)


def unregister():
    print("Template Unregister", __file__)
    import bpy

    register = bpy.utils.register_class
    for cls in class_store:
        register(cls)
    class_store.clear()

