
import sys
import platform


def parse_requirements(filename, platform, machine):
    from pip_requirements_parser import RequirementsFile
    import pkg_resources

    reqs = RequirementsFile.from_file(filename).to_dict()["requirements"]

    # monkey-patch in the platform and machine
    pkg_resources.sys.platform = platform

    def m():
        return machine

    pkg_resources.platform.machine = m

    evaluate_marker = pkg_resources.evaluate_marker

    deps = {}

    for req in reqs:
        name = req["name"]
        specifier = req["specifier"]
        marker = req["marker"]

        if len(specifier) == 0:
            specifier = ""
        else:
            specifier = specifier[0]

        if marker is not None:
            # check to see if this line fits this platform
            include = evaluate_marker(marker)
        else:
            include = True

        if include:
            deps[name] = specifier

    reqs = list(deps.keys())
    reqs.sort()

    result = []

    for req in reqs:
        result.append(f"{req}{deps[req]}")

    return result


try:
    reqs = sys.argv[1]
except Exception:
    print("USAGE: python test_requirements.py requirements_file.txt")
    sys.exit(-1)

current_platform = sys.platform
current_machine = platform.machine()

print(f"Running on {current_platform} : {current_machine}")

# I've checked these platforms - there is no agreement
# on what to call the processor...!
platforms = [ ("darwin", "arm64"),
              ("darwin", "x86_64"),
              ("linux", "aarch64"),
              ("linux", "x86_64"),
              ("win32", "AMD64") ]

for p in platforms:
    print(f"\nTesting {p[0]} : {p[1]}")
    result = parse_requirements(reqs, p[0], p[1])
    print(result)
