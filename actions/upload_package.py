import os
import sys
import glob

script = os.path.abspath(sys.argv[0])

try:
    channel = sys.argv[1]
except:
    channel = "dev"

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"sire source is in {srcdir}\n")
print(f"upload channel is {channel}\n")

# Get the anaconda token to authorise uploads
if "ANACONDA_TOKEN" in os.environ:
    conda_token = os.environ["ANACONDA_TOKEN"]
else:
    conda_token = "TEST"

# rattler-build outputs to an 'output' directory by default.
# Try rattler-build output first, fall back to conda-bld.
output_dir = os.path.join(srcdir, "output")

sire_pkg = glob.glob(os.path.join(output_dir, "**", "sire-*.conda"), recursive=True)

if len(sire_pkg) == 0:
    # Fall back to conda-bld location (legacy)
    if "CONDA" in os.environ:
        conda_bld = os.path.join(os.environ["CONDA"], "envs", "sire_build", "conda-bld")
        sire_pkg = glob.glob(os.path.join(conda_bld, "*-*", "sire-*.tar.bz2"))

    if len(sire_pkg) == 0:
        print("No sire packages to upload?")
        sys.exit(-1)

packages = sire_pkg

print(f"Uploading packages:")
print(" * ", "\n *  ".join(packages))

packages = " ".join(packages)


def run_cmd(cmd):
    import subprocess

    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()


gitdir = os.path.join(srcdir, ".git")

tag = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} tag --contains")

# The channel has been specified as an extra argument
# This will either be 'test' for main releases, or
# 'dev' for dev releases (the default)
channel = channel.lstrip().rstrip().replace(" ", "_").lower()

if len(channel) == 0:
    channel = "dev"

print(f"\nUploading to the '{channel}' channel")

label = f"--label {channel}"

# Upload the packages to the openbiosim channel on Anaconda Cloud.
cmd = f"anaconda --token {conda_token} upload --user openbiosim {label} --force {packages}"

print(f"\nUpload command:\n\n{cmd}\n")

# Label release packages with main and dev so that dev is at least as new as
# main.
if conda_token == "TEST":
    print("Not uploading as the ANACONDA_TOKEN is not set!")
    sys.exit(-1)

output = run_cmd(cmd)

print(output)

print("Package uploaded!")
