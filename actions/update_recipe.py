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
    print(f"Using environment from {sys.argv[1]}")

    env_channels = d["channels"]
else:
    env_reqs = []
    env_channels = []

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


def combine(reqs0, reqs1):
    """
    Combine requirements together, removing those from reqs0
    that appear in reqs1 (reqs1 has priority)
    """
    if type(reqs0) is not list:
        reqs0 = [reqs0]

    if type(reqs1) is not list:
        reqs1 = [reqs1]

    import re
    r = re.compile(r'([\w\d\-_]*)(>=|<=|==|<|>|=)(\d\.?\*?)*,?(>=|<=|=|==|<|>?)(\d\.?\*?)*|(\d\.?\*?)*\|(\d\.?\*?)*|(\d\.?\*?)*')

    reqs = []

    for req0 in reqs0:
        found = False

        m = r.match(req0)

        if m.groups()[0] is None:
            r0 = req0
        else:
            r0 = m.groups()[0]

        for req1 in reqs1:
            m = r.match(req1)

            if m.groups()[0] is None:
                req = req1
            else:
                req = m.groups()[0]

            if r0 == req:
                found = True
                break

        if not found:
            reqs.append(req0)

    return reqs + reqs1


def check_reqs(reqs0, reqs1):
    """
    Update reqs0 so that if there are any version requirements
    in reqs1 that affect dependencies in reqs0, then 
    reqs0 is updated to include those versions.
    """
    if type(reqs0) is not list:
        reqs0 = [reqs0]

    if type(reqs1) is not list:
        reqs1 = [reqs1]
    
    import re
    r = re.compile(r'([\w\d\-_]*)(>=|<=|==|<|>|=)(\d\.?\*?)*,?(>=|<=|=|==|<|>?)(\d\.?\*?)*|(\d\.?\*?)*\|(\d\.?\*?)*|(\d\.?\*?)*')

    reqs = []
                
    for req0 in reqs0:
        found = False
        found_req = None    
                
        m = r.match(req0)
    
        if m.groups()[0] is None:
            r0 = req0
        else:
            r0 = m.groups()[0]

        for req1 in reqs1:
            m = r.match(req1)

            if m.groups()[0] is None:
                req = req1
            else:
                req = m.groups()[0]
            
            if r0 == req:
                found = True
                found_req = req1
                break
            
        if found:
            reqs.append(found_req)
        else:
            reqs.append(req0)
        
    return reqs


build_reqs = dep_lines(check_reqs(build_reqs, env_reqs))
host_reqs = dep_lines(combine(host_reqs, env_reqs))
run_reqs = dep_lines(combine(run_reqs + bss_reqs, env_reqs))
test_reqs = dep_lines(test_reqs)

with open(recipe, "w") as FILE:
    for line in lines:
        if line.find("SIRE_BUILD_REQUIREMENTS") != -1:
            line = build_reqs
        elif line.find("SIRE_HOST_REQUIREMENTS") != -1:
            line = host_reqs
        elif line.find("SIRE_RUN_REQUIREMENTS") != -1:
            line = run_reqs
        elif line.find("SIRE_TEST_REQUIREMENTS") != -1:
            line = test_reqs
        else:
            line = line.replace("SIRE_REMOTE", sire_remote)
            line = line.replace("SIRE_BRANCH", sire_branch)

        FILE.write(line)

channels = ["conda-forge", "openbiosim/label/dev"]

for channel in env_channels:
    if channel not in channels:
        channels.insert(0, channel)

channels = " ".join([f"-c {x}" for x in channels])

print("\nBuild this package using the command")
print(f"conda mambabuild {channels} {condadir}")



