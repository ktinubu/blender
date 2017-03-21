
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

import bpy
from bpy.types import (
    Panel,
)

from bl_ui.properties_render import RenderButtonsPanel


# Just an example of a render panel
class RENDER_PT_render_simple(Panel, RenderButtonsPanel):
    bl_label = "101 Render Panel"
    COMPAT_ENGINES = {'BLENDER_RENDER'}
    def draw(self, context):
        layout = self.layout

        row = layout.row(align=True)
        row.operator("render.render", text="Render", icon='RENDER_STILL')


classes = (
    RENDER_PT_render_simple,
)

def register():
    print("Template UI Register", __file__)
    from bpy.utils import register_class

    for cls in classes:
        register_class(cls)


def unregister():
    print("Template UI Unregister", __file__)
    from bpy.utils import unregister_class

    for cls in classes:
        unregister_class(cls)