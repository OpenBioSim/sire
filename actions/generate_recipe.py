"""Generate a rattler-build recipe.yaml from pixi.toml.

This script reads the pixi.toml file (the single source of truth for
dependencies) and generates a rattler-build recipe.yaml with the
appropriate if/then conditional blocks for platform-specific dependencies.

Usage:
    python actions/generate_recipe.py [--features obs emle]

The --features flag controls which optional dependency groups are
included in the host section of the recipe.
"""

import argparse
import os
import subprocess
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


# Categorisation of core dependencies into conda recipe sections.
# Dependencies listed here go into "build"; everything else from
# [dependencies] goes into "host". The "run" list is a subset of
# host that is also needed at runtime.
BUILD_DEPS = {
    "cmake",
    "git",
    "make",
    "libtool",
    "pybind11",
    "sysroot_linux-64",
}

RUN_DEPS = {
    "gsl",
    "lazy_import",
    "libnetcdf",
    "openmm",
    "pandas",
    "python",
    "qt-main",
    "rich",
    "tbb",
}

# Mapping from pixi platform strings to rattler-build selector expressions.
PLATFORM_SELECTORS = {
    "linux-64": "linux and x86_64",
    "linux-aarch64": "linux and aarch64",
    "osx-arm64": "osx and arm64",
    "win-64": "win",
}

ALL_PLATFORMS = list(PLATFORM_SELECTORS.keys())


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a rattler-build recipe.yaml from pixi.toml"
    )
    parser.add_argument(
        "--features",
        nargs="*",
        default=[],
        help="Optional feature groups to include (e.g. obs emle)",
    )
    parser.add_argument(
        "--pixi-toml",
        default=None,
        help="Path to pixi.toml (default: auto-detect from repo root)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output path for recipe.yaml (default: recipes/sire/recipe.yaml)",
    )
    return parser.parse_args()


def run_cmd(cmd):
    """Run a shell command and return stripped stdout."""
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()


def get_git_info(srcdir):
    """Get the git remote URL and branch/tag."""
    gitdir = os.path.join(srcdir, ".git")

    remote = run_cmd(
        f"git --git-dir={gitdir} --work-tree={srcdir} config --get remote.origin.url"
    )
    if not remote.endswith(".git"):
        remote += ".git"

    branch = run_cmd(
        f"git --git-dir={gitdir} --work-tree={srcdir} rev-parse --abbrev-ref HEAD"
    )

    if branch == "HEAD":
        # Handle detached HEAD (e.g. in GitHub Actions PR checkouts).
        # GITHUB_HEAD_REF is set for pull_request events.
        # GITHUB_REF_NAME is set for push/workflow_dispatch events.
        branch = os.environ.get("GITHUB_HEAD_REF") or os.environ.get(
            "GITHUB_REF_NAME", ""
        )

        if not branch:
            # Fall back to tag detection for release builds.
            branch = run_cmd(
                f"git --git-dir={gitdir} --work-tree={srcdir} describe --tags"
            )
            if "-" in branch:
                raise RuntimeError("Cannot perform a tag build from a non-tag commit!")

    # Get the most recent tag for specifying the version in the recipe. If there are
    # no tags, use "PR".
    version = run_cmd(
        f"git --git-dir={gitdir} --work-tree={srcdir} describe --tags --abbrev=0"
    )
    if not version:
        version = "PR"

    # Work out the build number as the number of commits since the most recent tag,
    # or "0" if there are no commits.
    build = run_cmd(
        f"git --git-dir={gitdir} --work-tree={srcdir} rev-list --count {version}.."
    )
    if not build:
        build = "0"

    return remote, branch, version, build


def load_pixi_toml(path):
    """Load and parse the pixi.toml file."""
    with open(path, "rb") as f:
        return tomllib.load(f)


def format_dep(name, spec):
    """Format a dependency as a conda-style string."""
    if spec == "*" or spec == "":
        return name
    # Handle version specs that already start with an operator
    if spec[0] in ">=<!":
        return f"{name} {spec}"
    return f"{name} {spec}"


def extract_deps(data, section_path):
    """Extract dependencies from a nested dict path like 'target.linux-64.dependencies'.

    Returns a dict of {name: spec} or empty dict if the path doesn't exist.
    """
    parts = section_path.split(".")
    d = data
    for part in parts:
        if isinstance(d, dict) and part in d:
            d = d[part]
        else:
            return {}
    if not isinstance(d, dict):
        return {}
    return d


def categorise_deps(deps):
    """Split dependencies into build, host, and run categories.

    Returns (build_deps, host_deps, run_deps) as lists of formatted strings.
    """
    build = []
    host = []
    run = []

    for name, spec in sorted(deps.items()):
        formatted = format_dep(name, spec)
        if name in BUILD_DEPS:
            build.append(formatted)
        else:
            host.append(formatted)
            if name in RUN_DEPS:
                run.append(formatted)

    return build, host, run


def get_platform_deps(data, feature=None):
    """Get platform-specific deps, returning {platform: {name: spec}}.

    If feature is None, reads from top-level target sections.
    If feature is a string, reads from feature.{feature}.target sections.
    """
    result = {}
    for platform in ALL_PLATFORMS:
        if feature:
            path = f"feature.{feature}.target.{platform}.dependencies"
        else:
            path = f"target.{platform}.dependencies"
        result[platform] = extract_deps(data, path)
    return result


def deps_to_yaml_lines(deps, indent=4):
    """Convert a list of dep strings to indented YAML list lines."""
    prefix = " " * indent
    return [f"{prefix}- {dep}" for dep in sorted(deps)]


def platform_conditional_block(platform, deps, indent=4):
    """Generate an if/then block for platform-specific dependencies."""
    if not deps:
        return []

    prefix = " " * indent
    selector = PLATFORM_SELECTORS[platform]
    lines = [f"{prefix}- if: {selector}"]
    lines.append(f"{prefix}  then:")
    for dep in sorted(deps):
        lines.append(f"{prefix}    - {dep}")
    return lines


def generate_recipe(data, features, git_remote, git_branch, git_version, git_number):
    """Generate the complete recipe.yaml content."""
    # Extract core (all-platform) dependencies
    core_deps = extract_deps(data, "dependencies")

    # Extract platform-specific core deps
    platform_core = get_platform_deps(data)

    # Merge platform-specific into core for categorisation (all-platform deps)
    build_deps, host_deps, run_deps = categorise_deps(core_deps)

    # Collect platform-specific build and host deps
    platform_build = {}  # {platform: [dep_strings]}
    platform_host = {}
    for platform in ALL_PLATFORMS:
        pdeps = platform_core[platform]
        pb, ph, _ = categorise_deps(pdeps)
        platform_build[platform] = pb
        platform_host[platform] = ph

    # Collect feature deps (all go to host)
    feature_all_platform = []  # all-platform feature deps
    feature_platform = {p: [] for p in ALL_PLATFORMS}  # platform-specific feature deps

    for feature in features:
        # All-platform feature deps
        fdeps = extract_deps(data, f"feature.{feature}.dependencies")
        for name, spec in sorted(fdeps.items()):
            formatted = format_dep(name, spec)
            if formatted not in feature_all_platform:
                feature_all_platform.append(formatted)

        # Platform-specific feature deps
        fp = get_platform_deps(data, feature=feature)
        for platform in ALL_PLATFORMS:
            for name, spec in sorted(fp[platform].items()):
                formatted = format_dep(name, spec)
                if formatted not in feature_platform[platform]:
                    feature_platform[platform].append(formatted)

    # Also add feature all-platform deps to host
    for dep in feature_all_platform:
        if dep not in host_deps:
            host_deps.append(dep)

    # Extract test deps
    test_deps = []
    test_platform = {p: [] for p in ALL_PLATFORMS}
    tdeps = extract_deps(data, "feature.test.dependencies")
    for name, spec in sorted(tdeps.items()):
        test_deps.append(format_dep(name, spec))

    tp = get_platform_deps(data, feature="test")
    for platform in ALL_PLATFORMS:
        for name, spec in sorted(tp[platform].items()):
            test_platform[platform].append(format_dep(name, spec))

    # Build the YAML
    lines = []

    # Context
    lines.append("context:")
    lines.append("  name: sire")
    lines.append("")

    # Package
    lines.append("package:")
    lines.append("  name: ${{ name }}")
    lines.append(f"  version: {git_version}")
    lines.append("")

    # Source
    lines.append("source:")
    lines.append(f"  git: {git_remote}")
    lines.append(f"  branch: {git_branch}")
    lines.append("")

    # Build
    lines.append("build:")
    lines.append(f"  number: {git_number}")
    lines.append("  files:")
    lines.append("    exclude:")
    lines.append("      - etc/OpenCL/")
    lines.append("")

    # Requirements
    lines.append("requirements:")

    # Build requirements
    lines.append("  build:")
    lines.append("    - ${{ compiler('c') }}")
    lines.append("    - ${{ compiler('cxx') }}")
    lines.extend(deps_to_yaml_lines(build_deps, indent=4))
    for platform in ALL_PLATFORMS:
        lines.extend(
            platform_conditional_block(platform, platform_build[platform], indent=4)
        )

    # Host requirements
    lines.append("  host:")
    lines.extend(deps_to_yaml_lines(sorted(set(host_deps)), indent=4))
    # Platform-specific host deps (core + feature)
    for platform in ALL_PLATFORMS:
        combined = sorted(set(platform_host[platform] + feature_platform[platform]))
        lines.extend(platform_conditional_block(platform, combined, indent=4))

    # Run requirements
    lines.append("  run:")
    lines.extend(deps_to_yaml_lines(run_deps, indent=4))

    # Run constraints
    lines.append("  run_constraints:")
    lines.append("    - ${{ pin_compatible('gemmi', upper_bound='x.x.x') }}")
    lines.append("    - ${{ pin_compatible('openmm', upper_bound='x.x') }}")
    lines.append("    - ${{ pin_compatible('rdkit', upper_bound='x.x.x') }}")

    lines.append("")

    # Tests
    lines.append("tests:")

    # Import test
    lines.append("  - python:")
    lines.append("      imports:")
    for module in [
        "sire",
        "sire.analysis",
        "sire.base",
        "sire.cas",
        "sire.cluster",
        "sire.error",
        "sire.ff",
        "sire.id",
        "sire.io",
        "sire.maths",
        "sire.mm",
        "sire.mol",
        "sire.move",
        "sire.qt",
        "sire.squire",
        "sire.stream",
        "sire.system",
        "sire.units",
        "sire.vol",
    ]:
        lines.append(f"        - {module}")

    # Script test (pytest)
    lines.append("  - script:")
    lines.append("      - pytest -vvv --color=yes --runveryslow ./tests")
    lines.append("    files:")
    lines.append("      source:")
    lines.append("        - tests/")
    lines.append("    requirements:")
    lines.append("      run:")
    lines.extend(deps_to_yaml_lines(test_deps, indent=8))
    for platform in ALL_PLATFORMS:
        lines.extend(
            platform_conditional_block(platform, test_platform[platform], indent=8)
        )

    lines.append("")

    # About
    lines.append("about:")
    lines.append("  homepage: https://github.com/openbiosim/sire")
    lines.append("  license: GPL-3.0-or-later")
    lines.append("  license_file: LICENSE")
    lines.append('  summary: "An advanced molecular modelling framework."')
    lines.append("  repository: https://github.com/openbiosim/sire")
    lines.append("  documentation: https://sire.openbiosim.org")
    lines.append("")

    # Extra
    lines.append("extra:")
    lines.append("  recipe-maintainers:")
    lines.append("    - chryswoods")
    lines.append("    - jmichel80")
    lines.append("    - lohedges")

    return "\n".join(lines) + "\n"


def main():
    args = parse_args()

    # Determine paths
    script = os.path.abspath(sys.argv[0])
    srcdir = os.path.dirname(os.path.dirname(script))

    pixi_path = args.pixi_toml or os.path.join(srcdir, "pixi.toml")
    output_path = args.output or os.path.join(srcdir, "recipes", "sire", "recipe.yaml")

    print(f"Reading dependencies from {pixi_path}")
    print(f"Features: {args.features or '(none)'}")

    data = load_pixi_toml(pixi_path)
    git_remote, git_branch, git_version, git_number = get_git_info(srcdir)

    print(f"Git remote: {git_remote}")
    print(f"Git branch: {git_branch}")
    print(f"Git version: {git_version}")
    print(f"Git build number: {git_number}")

    recipe = generate_recipe(
        data, args.features, git_remote, git_branch, git_version, git_number
    )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(recipe)

    print(f"\nRecipe written to {output_path}")
    print("\nBuild this package using the command:")
    print(
        f"  rattler-build build --recipe {os.path.dirname(output_path)} "
        f"-c conda-forge -c openbiosim/label/dev"
    )


if __name__ == "__main__":
    main()
