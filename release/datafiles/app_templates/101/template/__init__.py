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
        # setup_classes
        "class_store",
        # setup_addons
        "sys_path",
        # setup_ui_filter
        "ui_filter_store"
    )

    _template_addons = (
        "template_101_library",
    )

    def __init__(self):
        self.class_store = []
        self.sys_path = []
        self.ui_filter_store = []

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

    def setup_ui_filter(self):
        import bl_app_override

        def filter_operator(op_id):
            return op_id not in {
                "transform.mirror",
                "sound.mixdown",
                "object.modifier_add",
                "object.forcefield_toggle",
            }

        def filter_property(ty, prop):
            return (ty, prop) not in {
                ("Object", "location"),
                ("Object", "rotation_euler"),
                ("Object", "scale"),
                ("RenderSettings", "filter_size"),
                ("RenderSettings", "frame_map_new"),
                ("RenderSettings", "frame_map_old"),
                ("RenderSettings", "pixel_aspect_x"),
                ("RenderSettings", "pixel_aspect_y"),
                ("RenderSettings", "pixel_filter_type"),
                ("RenderSettings", "use_border"),
                ("RenderSettings", "use_crop_to_border"),
                ("RenderSettings", "use_placeholder"),
                ("RenderSettings", "use_render_cache"),
                ("Scene", "frame_step"),
            }

        def filter_label(text):
            # print(text)
            return text not in {
                "Aspect Ratio:",
                "Time Remapping:",
            }

        self.ui_filter_store = bl_app_override.ui_draw_filter_register(
            filter_operator=filter_operator,
            filter_property=filter_property,
            filter_label=filter_label,
        )

    def teardown_ui_filter(self):
        import bl_app_override
        bl_app_override.ui_draw_filter_unregister(
            self.ui_filter_store
        )
        self.ui_filter_store = None


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
    app_state.setup_ui_filter()

    ui.register()

def unregister():
    print("Template Unregister", __file__)

    ui.unregister()

    app_state.teardown_classes()
    app_state.teardown_addons()
    app_state.teardown_ui_filter()
