
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
    white_list = kw.pop("white_list", None)
    black_list = kw.pop("black_list", None)
    kw_items = tuple(kw.items())
    for cls in cls_parent.__subclasses__():
        # same as is_registered()
        if "bl_rna" in cls.__dict__:
            if black_list is not None and cls.__name__ in black_list:
                continue
            if ((white_list is not None and cls.__name__ is white_list) or
                    all((getattr(cls, attr) in expect) for attr, expect in kw_items)):
                yield cls


def ui_draw_filter_register(
    *,
    classes=None,
    filter_operator=None,
    filter_property=None,
    filter_label=None,
):
    import bpy

    UILayout = bpy.types.UILayout

    if classes is None:
        classes = (
            bpy.types.Panel,
            bpy.types.Menu,
            bpy.types.Header,
        )

    class OperatorProperties_Fake:
        pass

    class UILayout_Fake(bpy.types.UILayout):
        __slots__ = ()

        def __getattribute__(self, attr):
            # ensure we always pass down UILayout_Fake instances
            if attr in {"row", "split", "column", "box", "column_flow"}:
                real_func = UILayout.__getattribute__(self, attr)

                def dummy_func(*args, **kw):
                    # print("wrapped", attr)
                    ret = real_func(*args, **kw)
                    return UILayout_Fake(ret)
                return dummy_func

            elif attr in {"operator", "operator_menu_enum", "operator_enum"}:
                if filter_operator is None:
                    return UILayout.__getattribute__(self, attr)

                real_func = UILayout.__getattribute__(self, attr)

                def dummy_func(*args, **kw):
                    # print("wrapped", attr)
                    if filter_operator(args[0]):
                        ret = real_func(*args, **kw)
                    else:
                        # UILayout.__getattribute__(self, "label")()
                        # may need to be set
                        ret = OperatorProperties_Fake()
                    return ret
                return dummy_func
            elif attr in {"prop", "prop_enum"}:
                if filter_property is None:
                    return UILayout.__getattribute__(self, attr)

                real_func = UILayout.__getattribute__(self, attr)

                def dummy_func(*args, **kw):
                    # print("wrapped", attr)
                    if filter_property(args[0].__class__.__name__, args[1]):
                        ret = real_func(*args, **kw)
                    else:
                        ret = None
                    return ret
                return dummy_func
            elif attr == "label":
                if filter_label is None:
                    return UILayout.__getattribute__(self, attr)

                real_func = UILayout.__getattribute__(self, attr)

                def dummy_func(*args, **kw):
                    # print("wrapped", attr)
                    if filter_label(args[0] if args else kw["text"]):
                        ret = real_func(*args, **kw)
                    else:
                        # ret = real_func()
                        ret = None
                    return ret
                return dummy_func
            else:
                return UILayout.__getattribute__(self, attr)
            # print(self, attr)

        def operator(*args, **kw):
            return super().operator(*args, **kw)

        def label(*args, **kw):
            text = args[1] if args else kw["text"]
            if filter_label(text):
                return super().label(*args, **kw)
            else:
                return super().label(args[0], "")


    def draw_override(func_orig, self_real, context):
        # simple, no wrapping
        # return func_orig(self_wrap, context)

        class Wrapper(self_real.__class__):
            __slots__ = ()
            def __getattribute__(self, attr):
                if attr == "layout":
                    return UILayout_Fake(self_real.layout)
                else:
                    cls = super()
                    try:
                        return cls.__getattr__(self, attr)
                    except AttributeError:
                        # class variable
                        return getattr(cls, attr)

            @property
            def layout(self):
                # print("wrapped")
                return self_real.layout

        return func_orig(Wrapper(self_real), context)

    ui_filter_store = []

    for cls in classes:
        for subcls in list(cls.__subclasses__()):
            if "draw" in subcls.__dict__:  # don't want to get parents draw()

                def replace_draw():
                    # function also serves to hold draw_old in a local name-space
                    draw_orig = subcls.draw

                    def draw(self, context):
                        return draw_override(draw_orig, self, context)
                    subcls.draw = draw

                ui_filter_store.append((subcls, "draw", subcls.draw))

                replace_draw()

    return ui_filter_store


def ui_draw_filter_unregister(ui_filter_store):
    for (obj, attr, value) in ui_filter_store:
        setattr(obj, attr, value)
