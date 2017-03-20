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

bl_info = {
    "name": "Blender 101 Library",
    "author": "None",
    "version": (1, 0),
    "blender": (2, 75, 0),
    "location": "View3D > Add > Mesh > New Object",
    "description": "Simple test case library",
    "warning": "",
    "wiki_url": "",
    "category": "Add Mesh",
}

import os

BASE_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")

import bpy
from bpy.types import Operator
from bpy.props import (
    FloatVectorProperty,
    StringProperty,
)
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector


class OBJECT_OT_add_object_from_blend(Operator, AddObjectHelper):
    """Create a new Mesh Object"""
    bl_idname = "mesh.add_object"
    bl_label = "Add Mesh Object"
    bl_options = {'REGISTER', 'UNDO'}

    scale = FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    filepath = StringProperty(
        name="Filepath",
        maxlen=1024,
        options={'HIDDEN'},
    )
    type_attr = StringProperty(
        name="Type",
        maxlen=64,
        options={'HIDDEN'},
    )

    def execute(self, context):

        name = self.name
        filepath = self.filepath
        type_attr = self.type_attr

        with bpy.data.libraries.load(filepath) as (data_from, data_to):
            setattr(data_to, type_attr, getattr(data_from, type_attr))

        for data in getattr(data_to, type_attr):
            object_data_add(context, data, operator=self)

        return {'FINISHED'}


def paths_from_id(dir_id):
    """(display_name, full_path)"""
    d = os.path.join(DATA_DIR, dir_id)
    for f in os.listdir(d):
        if f.endswith(".blend"):
            yield (bpy.path.display_name(f), os.path.join(d, f))


def add_object_menu(self, context):
    layout = self.layout
    for name, f in paths_from_id("object"):
        props = layout.operator(
            OBJECT_OT_add_object_from_blend.bl_idname,
            text=name,
            icon='PLUGIN',
        )
        props.filepath = f
        props.type_attr = "meshes"


def register():
    bpy.utils.register_class(OBJECT_OT_add_object_from_blend)
    bpy.types.INFO_MT_mesh_add.append(add_object_menu)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_object_from_blend)
    bpy.types.INFO_MT_mesh_add.remove(add_object_menu)


if __name__ == "__main__":
    register()
