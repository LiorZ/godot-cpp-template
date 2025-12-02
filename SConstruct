#!/usr/bin/env python
import os
import sys

from methods import print_error
from SCons.Script import ARGUMENTS


libname = "godotmol"
projectdir = "demo"

# Accept friendly target aliases.
_user_target = ARGUMENTS.get("target", "")
if _user_target == "release":
    ARGUMENTS["target"] = "template_release"
elif _user_target == "debug":
    ARGUMENTS["target"] = "template_debug"

localEnv = Environment(tools=["default"], PLATFORM="")

# Build profiles can be used to decrease compile times.
# You can either specify "disabled_classes", OR
# explicitly specify "enabled_classes" which disables all other classes.
# Modify the example file as needed and uncomment the line below or
# manually specify the build_profile parameter when running SCons.

# localEnv["build_profile"] = "build_profile.json"

customs = ["custom.py"]
customs = [os.path.abspath(path) for path in customs]

opts = Variables(customs, ARGUMENTS)
opts.Update(localEnv)

Help(opts.GenerateHelpText(localEnv))

env = localEnv.Clone()

if not (os.path.isdir("godot-cpp") and os.listdir("godot-cpp")):
    print_error("""godot-cpp is not available within this folder, as Git submodules haven't been initialized.
Run the following command to download godot-cpp:

    git submodule update --init --recursive""")
    sys.exit(1)

# Gemmi needs C++ exceptions enabled.
env["disable_exceptions"] = False

env = SConscript("godot-cpp/SConstruct", {"env": env, "customs": customs})

env.Append(CPPPATH=["src/"])
# Ensure exception support for Gemmi (godot-cpp disables it by default).
env.Append(CXXFLAGS=["-fexceptions"])

gemmi_include_dir = "gemmi/include"
if os.path.isdir(gemmi_include_dir):
    env.Append(CPPPATH=[gemmi_include_dir])
else:
    print("Warning: gemmi/include not found. Gemmi headers are required to build the extension.")

# Build and statically link Gemmi so the Godot extension doesn't depend on a
# system-wide Gemmi shared library at runtime.
gemmi_src_dir = "gemmi/src"
gemmi_sources = [
    f"{gemmi_src_dir}/align.cpp",
    f"{gemmi_src_dir}/assembly.cpp",
    f"{gemmi_src_dir}/calculate.cpp",
    f"{gemmi_src_dir}/ccp4.cpp",
    f"{gemmi_src_dir}/crd.cpp",
    f"{gemmi_src_dir}/ddl.cpp",
    f"{gemmi_src_dir}/eig3.cpp",
    f"{gemmi_src_dir}/fprime.cpp",
    f"{gemmi_src_dir}/gz.cpp",
    f"{gemmi_src_dir}/intensit.cpp",
    f"{gemmi_src_dir}/json.cpp",
    f"{gemmi_src_dir}/mmcif.cpp",
    f"{gemmi_src_dir}/mmread_gz.cpp",
    f"{gemmi_src_dir}/monlib.cpp",
    f"{gemmi_src_dir}/mtz.cpp",
    f"{gemmi_src_dir}/mtz2cif.cpp",
    f"{gemmi_src_dir}/pdb.cpp",
    f"{gemmi_src_dir}/polyheur.cpp",
    f"{gemmi_src_dir}/read_cif.cpp",
    f"{gemmi_src_dir}/resinfo.cpp",
    f"{gemmi_src_dir}/riding_h.cpp",
    f"{gemmi_src_dir}/select.cpp",
    f"{gemmi_src_dir}/sprintf.cpp",
    f"{gemmi_src_dir}/dssp.cpp",
    f"{gemmi_src_dir}/symmetry.cpp",
    f"{gemmi_src_dir}/to_json.cpp",
    f"{gemmi_src_dir}/to_mmcif.cpp",
    f"{gemmi_src_dir}/to_pdb.cpp",
    f"{gemmi_src_dir}/topo.cpp",
    f"{gemmi_src_dir}/xds_ascii.cpp",
]

gemmi_zlib_sources = [
    "gemmi/third_party/zlib/adler32.c",
    "gemmi/third_party/zlib/crc32.c",
    "gemmi/third_party/zlib/gzlib.c",
    "gemmi/third_party/zlib/gzread.c",
    "gemmi/third_party/zlib/inflate.c",
    "gemmi/third_party/zlib/inffast.c",
    "gemmi/third_party/zlib/inftrees.c",
    "gemmi/third_party/zlib/zutil.c",
]

gemmi_env = env.Clone()
gemmi_env.Append(
    CPPPATH=["gemmi/third_party"],
    CPPDEFINES=[
        "GEMMI_BUILD",
        "NO_GZCOMPRESS=1",
        "DYNAMIC_CRC_TABLE=1",
        "Z_HAVE_UNISTD_H=1",
    ],
    CCFLAGS=["-fPIC"],
)
gemmi_lib = gemmi_env.StaticLibrary("bin/{}/gemmi".format(env["platform"]), gemmi_sources + gemmi_zlib_sources)
env.Append(LIBS=[gemmi_lib])

# Gemmi relies on standard C++ exceptions.
env["use_exceptions"] = True
sources = Glob("src/*.cpp")

if env["target"] in ["editor", "template_debug"]:
    try:
        doc_data = env.GodotCPPDocData("src/gen/doc_data.gen.cpp", source=Glob("doc_classes/*.xml"))
        sources.append(doc_data)
    except AttributeError:
        print("Not including class reference as we're targeting a pre-4.3 baseline.")

# Normalize suffix to drop template markers.
# .dev doesn't inhibit compatibility, so we don't need to key it.
# .universal just means "compatible with all relevant arches" so we don't need to key it.
suffix = env["suffix"].replace(".dev", "").replace(".universal", "")
suffix = suffix.replace(".template_debug", ".debug").replace(".template_release", ".release")

lib_filename = "{}{}{}{}".format(env.subst('$SHLIBPREFIX'), libname, suffix, env.subst('$SHLIBSUFFIX'))

library = env.SharedLibrary(
    "bin/{}/{}".format(env['platform'], lib_filename),
    source=sources,
)

copy = env.Install("{}/bin/{}/".format(projectdir, env["platform"]), library)

default_args = [library, copy]
Default(*default_args)
