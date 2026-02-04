"""
Installation script for Sire

This assumes that the python that is used to execute this script
is within a conda / pixi environment that already has all of the
required dependencies installed. Use pixi to create the environment
before running this script, e.g.:

    pixi install -e dev

USAGE:

    python setup.py build              : Compile sire (takes a long time!)

    python setup.py install            : Will build sire and then install

    python setup.py install_module     : Will only install the Python module
"""

import glob
import json
import os
import platform
import subprocess
import shutil
import sys

try:
    # We have to check the version, but we can't do this by
    # importing sire (this will block installing of files on
    # windows, as we can't copy over a library file that is linked
    # into the application)
    curver = (
        subprocess.check_output(
            [sys.executable, os.path.join("build/get_sire_version.py")]
        )
        .decode("utf8")
        .lstrip()
        .rstrip()
    )

    if len(curver) == 0:
        curver = None
except Exception:
    curver = None

# If Sire is already installed, see if the version number has changed.
if curver:
    with open("version.txt", "r") as f:
        ver = f.readline().strip()

    # Abort if the versions differ.
    if curver != ver:
        raise EnvironmentError(
            f"This environment already contains an install of Sire version {curver}. "
            "Please delete the installation or create a new environment before "
            f"installing version {ver}. Also remove old build directories from 'build'."
        )

try:
    # Find out how much memory we have in total.
    # The wrappers need 4 GB per core to compile
    import psutil

    total_memory_gb = psutil.virtual_memory()[0] / (1024 * 1024 * 1024)
except Exception:
    total_memory_gb = None

# We can only run this script from the sire directory
curdir = os.path.abspath(".")

if os.path.abspath(os.path.dirname(sys.argv[0])) != curdir:
    print("You can only run this script from the sire directory")
    sys.exit(-1)

# Detect the environment prefix (conda or pixi)
if "PREFIX" in os.environ and "BUILD_PREFIX" in os.environ:
    print("This a build initiated by conda-build or rattler-build")
    conda_base = os.path.abspath(os.environ["PREFIX"])
    print(f"Setting conda-base to {conda_base}")
else:
    # Find the prefix from the Python executable
    conda_base = os.path.abspath(os.path.dirname(sys.executable))

    if os.path.basename(conda_base) == "bin":
        conda_base = os.path.dirname(conda_base)

if "CONDA_DEFAULT_ENV" in os.environ:
    conda_env = os.environ["CONDA_DEFAULT_ENV"]
else:
    conda_env = None

python_exe = None

if os.path.exists(os.path.join(conda_base, "python.exe")):
    # Windows
    conda_bin = os.path.join(conda_base, "Library", "bin")
    python_exe = os.path.join(conda_base, "python.exe")
elif os.path.exists(os.path.join(conda_base, "bin", "python")):
    # MacOS and Linux
    conda_bin = os.path.join(conda_base, "bin")
    python_exe = os.path.join(conda_bin, "python")
else:
    print(
        "Cannot find a 'python' binary in directory '%s'. "
        "Are you running this script using the python executable "
        "from a valid conda or pixi environment?" % conda_base
    )
    sys.exit(-1)


# Get the build operating system and processor
is_linux = False
is_windows = False
is_macos = False
platform_name = platform.system()
machine = platform.machine()

if machine in ["aarch64", "arm64"]:
    machine = "arm64"
elif machine in ["x86_64", "AMD64"]:
    machine = "x86_64"
else:
    print(f"Unrecognised architecture ({machine}). Compile at your own risk.")

exe_suffix = ""
if platform_name == "Linux":
    platform_name = "linux"
    is_linux = True
elif platform_name == "Darwin":
    platform_name = "macos"
    is_macos = True
elif platform_name == "Windows":
    platform_name = "windows"
    exe_suffix = ".exe"
    is_windows = True
else:
    print("Unrecognised build platform: %s" % platform_name)
    print("We cannot compile sire.")
    sys.exit(-1)

platform_string = f"{platform_name}_{machine}"


def parse_args():
    import argparse
    import multiprocessing

    parser = argparse.ArgumentParser()

    ncores = multiprocessing.cpu_count()

    if total_memory_gb is None:
        # it is safest to half the number of python build cores
        if ncores % 2 == 0:
            npycores = int(ncores / 2)
        else:
            npycores = int((ncores + 1) / 2)
    else:
        # we need at least 3 GB RAM per core
        npycores = min(ncores, int(total_memory_gb / 3))

        if npycores < 1:
            npycores = 1

        print(f"Total memory available: {total_memory_gb} GB")
        print(f"Number of python compile cores: {npycores}")

    parser.add_argument(
        "-C",
        "--corelib",
        action="append",
        nargs=1,
        metavar=("PARAMETER=VALUE",),
        default=[],
        help="pass CMake definitions for corelib",
    )
    parser.add_argument(
        "-W",
        "--wrapper",
        action="append",
        nargs=1,
        metavar=("PARAMETER=VALUE",),
        default=[],
        help="pass CMake definitions for wrapper",
    )
    parser.add_argument(
        "-G",
        "--generator",
        action="append",
        nargs=1,
        metavar=("GENERATOR",),
        default=[],
        help="pass CMake generator",
    )
    parser.add_argument(
        "-A",
        "--architecture",
        action="append",
        nargs=1,
        metavar=("ARCHITECTURE",),
        default=[],
        help="pass CMake generator architecture, e.g. WIN64",
    )
    parser.add_argument(
        "-n",
        "--ncores",
        action="store",
        type=int,
        nargs=1,
        metavar=("N_CORES",),
        default=[ncores],
        help="Number of CPU cores used for compiling corelib "
        "(defaults to all available on the current system)",
    )
    parser.add_argument(
        "-N",
        "--npycores",
        action="store",
        type=int,
        nargs=1,
        metavar=("N_PYTHON_CORES",),
        default=[npycores],
        help="Number of CPU cores used for compiling Python wrappers "
        "(defaults to the number of CPU cores used for compiling corelib)",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        default=False,
        help="Skip the build of the C++ code (only use if you know that "
        "the C++ code is already built)",
    )
    parser.add_argument(
        "--install-metadata",
        action="store_true",
        default=False,
        help="Install package metadata. This is useful when you are building "
        "from source but still want to be able to query the installation using "
        "conda list sire.",
    )
    parser.add_argument(
        "action",
        nargs="*",
        help="Should be one of 'build', 'install' or 'install_module'.\n"
        "\n [build] : Compile and install corelib, and just compile the wrappers.\n"
        " [install] : 'build' plus install the wrappers and install the module.\n"
        " [install_module] : Just install the module (no compilation).",
    )
    return parser.parse_args()


def add_default_cmake_defs(cmake_defs, ncores):
    for a in (
        "ANACONDA_BASE=%s" % conda_base.replace("\\", "/"),
        "BUILD_NCORES=%s" % ncores,
    ):
        found = False
        for d in cmake_defs:
            if a in d[0]:
                found = True
                break
        if not found:
            cmake_defs.append([a])

    if is_macos:
        # don't compile with AVX as the resulting binaries won't
        # work on M1 macs
        cmake_defs.append(["SIRE_DISABLE_AVX=ON"])
        cmake_defs.append(["SIRE_DISABLE_AVX512F=ON"])


def make_cmd(ncores, install=False):
    if is_windows:
        action = "INSTALL" if install else "ALL_BUILD"
        make_args = "%s -- /m:%s /p:Configuration=Release /p:Platform=x64" % (
            action,
            ncores,
        )
    else:
        action = "install" if install else "all"
        make_args = "%s -- VERBOSE=1 -j %s" % (action, ncores)

    return make_args.split()


def _get_build_ext():
    if "CONDA_BUILD" in os.environ and os.environ["CONDA_BUILD"] == "1":
        return "conda_build"
    else:
        if conda_env is not None:
            ext = "_" + conda_env.replace(" ", "_").replace(".", "_")
        else:
            ext = ""

        return os.path.basename(conda_base.replace(" ", "_").replace(".", "_")) + ext


def _get_bin_dir():
    if "CONDA_BUILD" in os.environ and os.environ["CONDA_BUILD"] == "1":
        bindir = os.environ["BUILD_PREFIX"]

        if is_windows:
            return os.path.join(bindir, "Library", "bin")
        else:
            return os.path.join(bindir, "bin")
    else:
        return conda_bin


def build(ncores: int = 1, npycores: int = 1, coredefs=[], pydefs=[]):
    print("\nCompiling the C++ code")

    CC = None
    CXX = None

    cmake = "cmake"

    if is_windows:
        cmake = f"{cmake}.exe"

    bindir = _get_bin_dir()
    cmake = os.path.join(bindir, cmake)

    if "CONDA_BUILD" in os.environ and os.environ["CONDA_BUILD"] == "1":
        conda_build = True
    else:
        os.environ["CONDA_BUILD"] = "0"
        conda_build = False

    # get the compilers
    if conda_build:
        print("This is a conda/rattler build")

        if is_windows:
            # Windows: vcvars is activated, let CMake find the MSVC.
            CXX = None
            CC = None
        else:
            # Try to get compilers from environment
            CXX = os.environ.get("CXX")
            CC = os.environ.get("CC")

    elif is_macos:
        try:
            CXX = glob.glob(os.path.join(bindir, "clang++"))[0]
            CC = glob.glob(os.path.join(bindir, "clang"))[0]
        except Exception:
            print("Cannot find the conda clang++ binaries!")
            print("Please ensure your environment has clang and clangxx installed.")
            print("If using pixi, run: pixi install -e dev")
            sys.exit(-1)

    elif is_linux:
        try:
            CXX = glob.glob(os.path.join(bindir, "*-g++"))[0]
            CC = glob.glob(os.path.join(bindir, "*-gcc"))[0]
        except Exception:
            print("Cannot find the conda g++ binaries!")
            print("Please ensure your environment has gcc and gxx installed.")
            print("If using pixi, run: pixi install -e dev")
            sys.exit(-1)

    if CC is not None and CXX is not None:
        print(f"Using compilers {CC} | {CXX}")

    # Make sure all of the above output is printed to the screen
    # before we start running any actual compilation
    sys.stdout.flush()

    # Now that the dependencies are installed, the next step
    # is to use cmake to build the corelib and wrapper in the build/corelib
    # and build/wrapper directories

    # change into the build/corelib directory
    OLDPWD = os.getcwd()

    build_ext = _get_build_ext()

    build_dir = os.path.abspath("build")

    coredir = os.path.join(build_dir, f"{build_ext}_corelib")

    if not os.path.exists(coredir):
        os.makedirs(coredir)

    if not os.path.isdir(coredir):
        print("SOMETHING IS WRONG. %s is not a directory?" % coredir)
        sys.exit(-1)

    os.chdir(coredir)

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = subprocess.run([cmake, "."])
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(
            os.path.dirname(os.path.dirname(os.getcwd())), "corelib"
        )

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print(
                "SOMETHING IS WRONG. There is no file %s"
                % os.path.join(sourcedir, "CMakeLists.txt")
            )
            sys.exit(-1)

        for a in ("NetCDF_ROOT_DIR", "OPENMM_ROOT_DIR"):
            for i, d in enumerate(args.corelib):
                if a in d[0]:
                    v = args.corelib.pop(i)[0]
                    if not a in os.environ:
                        os.environ[a] = v.split("=")[-1]

        add_default_cmake_defs(args.corelib, ncores)

        cmake_cmd = [
            cmake,
            *sum([["-D", d[0]] for d in args.corelib], []),
            *sum([["-G", g[0]] for g in args.generator], []),
            *sum([["-A", a[0]] for a in args.architecture], []),
            sourcedir,
        ]

        if CC is not None:
            os.environ["CC"] = CC
        if CXX is not None:
            os.environ["CXX"] = CXX

        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON CORELIB!")
        print("\n== OUTPUT LOG ==")
        try:
            for line in open("CMakeFiles/CMakeOutput.log").readlines():
                print(line.strip())
        except Exception as e:
            print(e)

        print("\n== ERROR LOG ==")
        try:
            for line in open("CMakeFiles/CMakeError.log").readlines():
                print(line.strip())
        except Exception as e:
            print(e)

        sys.exit(-1)

    # Now that cmake has run, we can compile and install corelib
    ######
    ###### Compiling and installing corelib
    ######
    # Compile and install, as need to install to compile the wrappers
    make_args = make_cmd(ncores, True)

    print('NOW RUNNING "%s" --build . --target %s' % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING CORELIB!")
        sys.exit(-1)

    ######
    ###### Compiling wrapper
    ######

    # Make sure that the Python in conda is used
    pydefs.append([f"PYTHON_EXECUTABLE={python_exe}"])

    os.chdir(OLDPWD)

    wrapperdir = os.path.join(build_dir, f"{build_ext}_wrapper")

    if not os.path.exists(wrapperdir):
        os.makedirs(wrapperdir)

    if not os.path.isdir(wrapperdir):
        print("SOMETHING IS WRONG. %s is not a directory?" % wrapperdir)
        sys.exit(-1)

    os.chdir(wrapperdir)

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = subprocess.run([cmake, "."])
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(
            os.path.dirname(os.path.dirname(os.getcwd())), "wrapper"
        )

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print(
                "SOMETHING IS WRONG. There is no file %s"
                % os.path.join(sourcedir, "CMakeLists.txt")
            )
            sys.exit(-1)

        add_default_cmake_defs(args.wrapper, npycores)

        cmake_cmd = [
            cmake,
            *sum([["-D", d[0]] for d in args.wrapper], []),
            *sum([["-G", g[0]] for g in args.generator], []),
            *sum([["-A", a[0]] for a in args.architecture], []),
            sourcedir,
        ]

        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON WRAPPER!")
        sys.exit(-1)

    # Just compile the wrappers
    make_args = make_cmd(npycores, False)

    print('NOW RUNNING "%s" --build . --target %s' % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING THE WRAPPERS!")
        sys.exit(-1)

    os.chdir(OLDPWD)


def install_module(ncores: int = 1):
    print("\nInstalling the module")

    OLDPWD = os.getcwd()

    build_ext = _get_build_ext()

    build_dir = os.path.abspath("build")

    cmake = "cmake"

    if is_windows:
        cmake = f"{cmake}.exe"

    bindir = _get_bin_dir()
    cmake = os.path.join(bindir, cmake)

    moduledir = os.path.join(build_dir, f"{build_ext}_module")

    if not os.path.exists(moduledir):
        os.makedirs(moduledir)

    if not os.path.isdir(moduledir):
        print("SOMETHING IS WRONG. %s is not a directory?" % moduledir)
        sys.exit(-1)

    os.chdir(moduledir)

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = subprocess.run([cmake, "."])
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(
            os.path.dirname(os.path.dirname(os.getcwd())), "src", "sire"
        )

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print(
                "SOMETHING IS WRONG. There is no file %s"
                % os.path.join(sourcedir, "CMakeLists.txt")
            )
            sys.exit(-1)

        add_default_cmake_defs(args.wrapper, ncores)
        cmake_cmd = [
            cmake,
            *sum([["-D", d[0]] for d in args.wrapper], []),
            *sum([["-G", g[0]] for g in args.generator], []),
            *sum([["-A", a[0]] for a in args.architecture], []),
            sourcedir,
        ]
        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON MODULE!")
        sys.exit(-1)

    make_args = make_cmd(ncores, True)

    # Now that cmake has run, we can compile and install wrapper
    print('NOW RUNNING "%s" --build . --target %s' % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
        sys.exit(-1)


def install(ncores: int = 1, npycores: int = 1):
    print("\nInstalling sire")

    OLDPWD = os.getcwd()

    build_ext = _get_build_ext()

    build_dir = os.path.abspath("build")

    wrapperdir = os.path.join(build_dir, f"{build_ext}_wrapper")

    cmake = "cmake"

    if is_windows:
        cmake = f"{cmake}.exe"

    bindir = _get_bin_dir()
    cmake = os.path.join(bindir, cmake)

    if not os.path.exists(wrapperdir):
        os.makedirs(wrapperdir)

    if not os.path.isdir(wrapperdir):
        print("SOMETHING IS WRONG. %s is not a directory?" % wrapperdir)
        sys.exit(-1)

    os.chdir(wrapperdir)

    # Now install the wrappers
    make_args = make_cmd(npycores, True)

    print('NOW RUNNING "%s" --build . --target %s' % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING THE WRAPPERS!")
        sys.exit(-1)

    os.chdir(OLDPWD)

    install_module(ncores=ncores)


if __name__ == "__main__":
    OLDPWD = os.getcwd()

    args = parse_args()

    if len(args.action) != 1:
        print("Please use either 'build', 'install' or 'install_module'")
        sys.exit(-1)

    action = args.action[0]

    if is_windows and (args.generator is None or len(args.generator) == 0):
        args.generator = [["Visual Studio 17 2022"]]
        args.architecture = [["x64"]]
    elif is_macos:
        # fix compile bug when INSTALL_NAME_TOOL is not set
        if "INSTALL_NAME_TOOL" not in os.environ:
            os.environ["INSTALL_NAME_TOOL"] = "install_name_tool"

    if action == "install":
        if not args.skip_build:
            build(
                ncores=args.ncores[0],
                npycores=args.npycores[0],
                coredefs=args.corelib,
                pydefs=args.wrapper,
            )

        install(ncores=args.ncores[0], npycores=args.npycores[0])

    elif action == "build":
        build(
            ncores=args.ncores[0],
            npycores=args.npycores[0],
            coredefs=args.corelib,
            pydefs=args.wrapper,
        )

    elif action == "install_module":
        install_module(ncores=args.ncores[0])

    else:
        print(
            f"Unrecognised action '{action}'. Please use "
            "'build', 'install' or 'install_module'"
        )

    # Create minimist package metadata so that 'conda list sire' works.
    if args.install_metadata:
        os.chdir(OLDPWD)
        if "CONDA_PREFIX" in os.environ:
            metadata_dir = os.path.join(os.environ["CONDA_PREFIX"], "conda-meta")
            if os.path.exists(metadata_dir):
                # Get the Python version.
                pyver = f"py{sys.version_info.major}{sys.version_info.minor}"
                metadata = {
                    "name": "sire",
                    "version": open("version.txt").readline().strip(),
                    "build": pyver,
                    "build_number": 0,
                    "channel": "local",
                    "size": 0,
                    "license": "GPL-3.0-or-later",
                    "subdir": platform_string,
                }
                metadata_file = os.path.join(
                    metadata_dir, f"sire-{metadata['version']}-{pyver}.json"
                )
                with open(metadata_file, "w") as f:
                    json.dump(metadata, f, indent=2)
                print(f"Created conda package metadata file: {metadata_file}")
