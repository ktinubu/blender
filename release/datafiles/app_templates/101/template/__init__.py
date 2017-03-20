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
import bl_app_override


class AppStateStore:
    # Utility class to encapsulate application state, backup and restore.
    __slots__ = (
        "class_store",
        "sys_path",
    )

    _template_addons = (
    )

    def __init__(self):
        self.class_store = []
        self.sys_path = []

    def setup_classes(self):
        assert(len(self.class_store) == 0)

        # Classes

        self.class_store.extend(
            bl_app_override.class_filter(
                bpy.types.Panel,
                # match any of these values
                bl_region_type={'TOOLS', 'WINDOW'},
                bl_space_type={'VIEW_3D', 'PROPERTIES'},
                # keep basic panels
                black_list={
                    'OBJECT_PT_transform',
                    'VIEW3D_PT_tools_add_object',
                    'VIEW3D_PT_tools_meshedit',
                },
            ),
        )

        unregister = bpy.utils.unregister_class
        for cls in self.class_store:
            unregister(cls)

    def teardown_classes(self):
        assert(len(self.class_store) != 0)

        register = bpy.utils.register_class
        for cls in self.class_store:
            register(cls)
        self.class_store.clear()

    def setup_addons(self):
        import sys
        import os
        template_addons = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "addons"))
        if template_addons not in sys.path:
            sys.path.append(template_addons)
            self.sys_path.append(template_addons)

        import addon_utils
        for addon in self._template_addons:
            addon_utils.enable(addon)

    def teardown_addons(self):
        import sys
        for path in self.sys_path:
            # should always succeed, but if not its no problem
            try:
                sys.path.remove(path)
            except:
                pass
        self.sys_path.clear()

        import addon_utils
        for addon in self._template_addons:
            addon_utils.disable(addon)


app_state = AppStateStore()

from . import ui

def register():
    print("Template Register", __file__)
    app_state.setup_classes()
    app_state.setup_addons()

    ui.register()

def unregister():
    print("Template Unregister", __file__)

    ui.unregister()

    app_state.teardown_classes()
    app_state.teardown_addons()
