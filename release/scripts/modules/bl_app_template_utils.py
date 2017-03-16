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
Similar to ``addon_utils``, except we can only have one active at a time.

In most cases users of this module will simply call 'activate'.
"""

__all__ = (
    # "paths",
    # "modules",
    # "check",
    # "enable",
    # "disable",

    "activate",
    "import_from_path",
    "import_from_id",
    "reset",
    # "module_bl_info",
)

import bpy as _bpy

# Normally matches 'user_preferences.app_template_id',
# but loading new preferences will get us out of sync.
_app_template = {
    "id": "",
}

# instead of sys.modules
# note that we only ever have one template enabled at a time
# so it may not seem necessary to use this.
#
# However, templates may want to share between each-other,
# so any loaded modules are stored here?
#
# Note that the ID here is the app_template_id , not the modules __name__.
_modules = {}


# -----------------------------------------------------------------------------
# Helper Classes

# Avoids leaving sys.paths & modules in an unknown state.
class _IsolateImportHelper:

    __slots__ = ("path", "module", "module_name")

    def __init__(self, path, module_name):
        self.path = path
        self.module_name = module_name

    def __enter__(self):
        import sys
        self.module = sys.modules.pop(self.module_name, None)
        sys.path.insert(0, self.path)

    def __exit__(self, type, value, traceback):
        import sys
        if self.module is not None:
            sys.modules[self.module_name] = self.module
        else:
            sys.modules.pop(self.module_name, None)
        try:
            sys.path.remove(self.path)
        except Exception:
            pass


def _enable(template_id, *, handle_error=None):
    import os
    import sys
    from bpy_restrict_state import RestrictBlend

    if handle_error is None:
        def handle_error(ex):
            import traceback
            traceback.print_exc()

    # Split registering up into 3 steps so we can undo
    # if it fails par way through.

    # disable the context, using the context at all is
    # really bad while loading an template, don't do it!
    with RestrictBlend():

        # 1) try import
        try:
            mod = import_from_id(template_id)
            mod.__template_enabled__ = False
            _modules[template_id] = mod
        except Exception as ex:
            handle_error(ex)
            return None

        # 2) try register collected modules
        # removed, templates need to handle own registration now.

        # 3) try run the modules register function
        try:
            mod.register()
        except Exception as ex:
            print("Exception in module register(): %r" %
                  getattr(mod, "__file__", template_id))
            handle_error(ex)
            del _modules[template_id]
            return None

    # * OK loaded successfully! *
    mod.__template_enabled__ = True

    if _bpy.app.debug_python:
        print("\tapp_template_utils.enable", mod.__name__)

    return mod


def _disable(template_id, *, handle_error=None):
    """
    Disables a template by name.

    :arg template_id: The name of the template and module.
    :type template_id: string
    :arg handle_error: Called in the case of an error, taking an exception argument.
    :type handle_error: function
    """
    import sys

    if handle_error is None:
        def handle_error(ex):
            import traceback
            traceback.print_exc()

    mod = _modules.get(template_id)

    # possible this addon is from a previous session and didn't load a
    # module this time. So even if the module is not found, still disable
    # the addon in the user prefs.
    if mod and getattr(mod, "__template_enabled__", False) is not False:
        mod.__template_enabled__ = False

        try:
            mod.unregister()
        except Exception as ex:
            print("Exception in module unregister(): %r" %
                  getattr(mod, "__file__", template_id))
            handle_error(ex)
    else:
        print("addon_utils.disable: %s not %s." %
              (template_id, "disabled" if mod is None else "loaded"))

    if _bpy.app.debug_python:
        print("\tapp_template_utils.disable", template_id)

def import_from_path(path):
    """
    Imports 'startup' from a path.
    """

    module_name = "template"
    # loading packages without modifying sys.path is some dark-art.
    # for now just use regular import but don't use sys.modules for cache.
    with _IsolateImportHelper(path, module_name):
        return __import__(module_name)


def import_from_id(template_id):
    path = next(iter(_bpy.utils.app_template_paths(template_id)), None)
    if path is None:
        raise Exception("%r template not found!" % template_id)
    else:
        return import_from_path(path)


def activate(template_id=None):

    # Disable all addons, afterwards caller must reset.
    import addon_utils
    addon_utils.disable_all()

    template_id_prev = _app_template["id"]
    if template_id_prev:
        _disable(template_id_prev)

    mod = _enable(template_id) if template_id else None

    if mod is not None:
        _app_template["id"] = template_id


def reset(*, reload_scripts=False):
    """
    Sets default state.
    """
    template_id = _bpy.context.user_preferences.app_template
    if _bpy.app.debug_python:
        print("bl_app_template_utils.reset('%s')" % template_id)

    if reload_scripts and False:
        # TODO, seems correct but reload fails
        import importlib
        import os
        import sys
        _modules_new = {}
        for key, mod in _modules.items():
            # Will always be 'template' but just use convention of __name__ to be sure.
            module_name = mod.__name__
            with _IsolateImportHelper(os.path.dirname(mod.__file__), module_name):
                sys.modules[module_name] = mod
                _modules_new[key] = importlib.reload(mod)
                del sys.modules[module_name]
        _modules.clear()
        _modules.update(_modules_new)
        del _modules_new

    activate(template_id)
