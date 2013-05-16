import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

m = {}

def browse_tree(root, folder, ext_modules):
    folder_path = os.path.join(root, folder)
    for entry in os.listdir(folder_path):
        fullpath = os.path.join(root, folder, entry)
        if os.path.isfile(fullpath):
            if folder == "":
                continue
            filename, extension = os.path.splitext(entry)
            if extension == ".py" and filename != "__init__":
                module_name = folder.replace("/",".") + "." + os.path.splitext(entry)[0]
                ext_modules.append(Extension(module_name, [fullpath]))
                m[module_name] = fullpath
        elif os.path.isdir(fullpath):
            if entry in ["timing", "tests", "build", ".git", "cython"]:
                continue
            browse_tree(root, os.path.join(folder, entry), ext_modules)

def create_modules():
    ext_modules = []
    current_folder = os.path.realpath("..")
    browse_tree(current_folder, "", ext_modules)
    if len(ext_modules):
        setup(
            name = 'pyPbrt',
            cmdclass = {'build_ext': build_ext},
            ext_modules = ext_modules
            )

create_modules()
