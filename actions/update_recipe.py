import sys
import os
import subprocess

script = os.path.abspath(sys.argv[0])

# we want to import the 'get_requirements' package from this directory
sys.path.insert(0, os.path.dirname(script))

from parse_requirements import parse_requirements

# has the user supplied an environment.yml file?
if len(sys.argv) > 1:
    from pathlib import Path
    import yaml
    d = yaml.safe_load(Path(sys.argv[1]).read_text())
    env_reqs = [x for x in d["dependencies"] if type(x) is str]
else:
    env_reqs = []

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"Sire source is in {srcdir}")

condadir = os.path.join(srcdir, "recipes", "sire")

print(f"conda recipe in {condadir}")

# Store the name of the recipe and template YAML files.
recipe = os.path.join(condadir, "meta.yaml")
template = os.path.join(condadir, "template.yaml")

# Now parse all of the requirements
build_reqs = parse_requirements(os.path.join(srcdir, "requirements_build.txt"))
print(build_reqs)
host_reqs = parse_requirements(os.path.join(srcdir, "requirements_host.txt"))
print(host_reqs)
run_reqs = parse_requirements(os.path.join(srcdir, "requirements_run.txt"))
print(run_reqs)
bss_reqs = parse_requirements(os.path.join(srcdir, "requirements_bss.txt"))
print(bss_reqs)
test_reqs = parse_requirements(os.path.join(srcdir, "requirements_test.txt"))
print(test_reqs)
print(env_reqs)


def run_cmd(cmd):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()


gitdir = os.path.join(srcdir, ".git")

# Get the sire remote
sire_remote = run_cmd(
    f"git --git-dir={gitdir} --work-tree={srcdir} config --get remote.origin.url"
)
sire_remote += ".git"
print(sire_remote)

# Get the Sire branch.
sire_branch = run_cmd(
    f"git --git-dir={gitdir} --work-tree={srcdir} rev-parse --abbrev-ref HEAD"
)
print(sire_branch)

lines = open(template, "r").readlines()


def dep_lines(deps):
    lines = []

    for dep in deps:
        lines.append(f"    - {dep}\n")

    return "".join(lines)

build_reqs = dep_lines(build_reqs)
host_reqs = dep_lines(host_reqs + env_reqs)
run_reqs = dep_lines(run_reqs + env_reqs)
bss_reqs = dep_lines(bss_reqs)
test_reqs = dep_lines(test_reqs)

with open(recipe, "w") as FILE:
    for line in lines:
        if line.find("SIRE_BUILD_REQUIREMENTS") != -1:
            line = build_reqs
        elif line.find("SIRE_HOST_REQUIREMENTS") != -1:
            line = host_reqs
        elif line.find("SIRE_RUN_REQUIREMENTS") != -1:
            line = run_reqs
        elif line.find("SIRE_BSS_REQUIREMENTS") != -1:
            line = bss_reqs
        elif line.find("SIRE_TEST_REQUIREMENTS") != -1:
            line = test_reqs
        else:
            line = line.replace("SIRE_REMOTE", sire_remote)
            line = line.replace("SIRE_BRANCH", sire_branch)

        FILE.write(line)
