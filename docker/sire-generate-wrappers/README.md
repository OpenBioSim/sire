# Generating the Python wrappers

This is a container that is used to consistently generate
the Python wrappers for sire. The container has a working
version of Py++ in a miniconda that also has all of the header
files that are needed for sire.

## Running the container

The simplest way to start is to run the container that we've
made already. Do this by typing;

```
$ docker run -it openbiosim/sire-generate-wrappers:x86
```

on Intel/AMD (X86-64) systems, and

```
$ docker run -it openbiosim/sire-generate-wrappers:arm64
```

on Apple / ARM64 systems.

Assuming you have docker installed and working, this should
pull the container from docker hub, and then start a bash
prompt in that container. You should see something like
this;

```
Unable to find image 'openbiosim/sire-generate-wrappers:arm64' locally
Trying to pull repository docker.io/siremol/sire-generate-wrappers ...
latest: Pulling from docker.io/siremol/sire-generate-wrappers
a2abf6c4d29d: Pull complete
c256cb8a03f5: Pull complete
96470ebef4ad: Pull complete
fe843007773c: Pull complete
13f39c13126f: Pull complete
0470e7ed57ae: Pull complete
5363ab184b58: Pull complete
d719333d5423: Pull complete
29c0367973b3: Pull complete
4642f7c70ef4: Pull complete
b70e0e58d6c8: Pull complete
Digest: sha256:a8a4513655dd43cebac306f77ac85a3fd953802759029bc771fe51db058eea95
Status: Downloaded newer image for siremol/sire-generate-wrappers:latest
(base) root:~#
```

## Generating the wrappers

You generate the wrappers for your chosen branch of sire by typing;

```
(base) root:~# generate_wrappers --branch {BRANCH_NAME}
```

where `{BRANCH_NAME}` is the name of the branch you want to use.
E.g. if you want to generate the wrappers for `devel`, you
would type

```
(base) root:~# generate_wrappers --branch devel
```

This will run in parallel, but be aware that this can take
a long time!

## Checking the wrappers

The `generate_wrappers` command will check out your branch
to the folder `$HOME/sire`. Feel free to explore this
folder and run commands manually if there are any problems.

You can check which wrappers were changed by running
`git status` in any of the wrapper directories.

## Pushing the wrappers to GitHub

You can push your new wrappers back to GitHub by
typing the command;

```
(base) root:~# push_wrappers --name "{YOUR NAME}" --email "{YOUR EMAIL}"
```

where `${YOUR NAME}` is your real name, as on GitHub, and
`${YOUR EMAIL}` is your real email. Note that you should put
these in double quotes.

This will use `git config` to set those values, before
running `git add` in `$HOME/sire`, and then running
`git commit` and `git push`.

The `git push` command will ask you for your GitHub
username and password. Your new wrappers will
be pushed if these are entered and are correct.

Assuming everything has worked, then congratulations!
You have successfully generated and pushed your new
sire Python wrappers.

## Using the container against an external sire directory

You may wish to use this container to create wrappers for
a sire directory on your computer. To do this, you can mount
the sire directory into the container via

```
docker run -it -v /path/to/sire:/root/sire openbiosim/sire-generate-wrappers:x86
```

or

```
docker run -it -v /path/to/sire:/root/sire openbiosim/sire-generate-wrappers:arm64
```

This will make your sire directory (in `/path/to/sire`) available in the
container as `/root/sire`. You can now generate wrappers in the
container that will write directly to your real local sire directory
(and so can be compiled, tested and pushed to GitHub from there).

Once in the container, you need to run `unpack_headers` to unpack the
headers.

```
(base) root:~# unpack_headers
```

You can now change into the sire directory and generate wrappers as you
need, e.g.

```
(base) root:~# cd sire/wrapper
(base) root:~/sire/wrapper# python AutoGenerate/scanheaders.py ~/sire/corelib/src/libs/ .
(base) root:~/sire/wrapper# cd Search
(base) root:~/sire/wrapper/Search# python ../AutoGenerate/create_wrappers.py
```

would regenerate the wrappers for the `Search` module. Note that this
will be slower than working entirely in the container due to the slowness
of translating between the container's filesystem and the host's
filesystem. But it is useful for iterative development as it avoids
you needing to push and pull changes via git.

Be very careful doing this, as the git version in the container is not
going to be the same as on your local computer. Also be careful not
to run `generate_wrappers` as this will delete your local sire container.

## Creating the sire-generate-wrappers container

You only need to read here if you want to create the
sire-generate-wrappers container yourself.

First, you must have a working installation on sire
on your computer. This is needed so that we can
copy the header files from this installation into
the container. We will assume that sire is installed
in `$HOME/sire.app`.

Next, you must navigate to this directory, e.g.

```
cd $HOME/sire/docker/sire-generate-wrappers
```

(assuming sire was cloned to `$HOME/sire`)

Create the `includes.tar.bz2` file by running the
`create_includes_tarball` script via;

```
./create_includes_tarball --sire $HOME/sire.app
```

Next, create the container via

```
docker build -t siremol/sire-generate-wrappers .
```

Assuming this worked, you can run the container via
the same command as at the top, e.g.

```
docker run -it siremol/sire-generate-wrappers
```

and the follow the rest of the instructions.

## Removing all docker containers

Sometimes you will want to remove ALL docker containers
from your computer. You can do this via these two
commands;

```
docker ps -a -q | xargs docker rm
docker images -a -q | xargs docker rmi -f
```
