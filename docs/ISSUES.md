# Fixing your genepred installation

Here you can find some troubleshooting in case your installation didn't go as planned.

## Cloning errors

None documented yet.

## Dependency errors

A common issue may be that you did not manage to install correctly all the dependencies. In that case, the error is not in this pipeline but in the installation itself. Some common issues are reported here.

### Using `./configure`, `make`, and `make install`

Most tools are installed with a sequence of commands:

```
./configure
make
make install
```

The first command `./configure` runs a script contained in the tool main directory called `configure`, which generates the so-called `Makefile` based on the system specifics.

The second command `make` executes the instructions contained in the `Makefile`.

The third command `make install` projects the installation onto the whole system. This command is not mandatory, as long as you use the tool specifying its full path to this directory. If you want it in your `$PATH` however, unless you add this directory to the `$PATH`, you will need to install it system-wide with this command.

Note that, on computing clusters, `make install` requires root permissions which you normally don't have. What people usually do is to have their own `local/bin` directory in which they have writing permissions, and add that `local/bin` directory to their won `$PATH`.

More about this here:
https://www.baeldung.com/linux/change-install-dir-make-install

### Tool installed, `check_dependencies.py` doesn't see it

This could be fixed very quickly. First, make sure of the exact name your installed dependency has in the `$PATH`. For example, it could be that your **bedtools** distribution in your cluster is called `bedtools_2.29` instead of `bedtools`.

Open the `check_dependencies.py` script:

```
nano check_dependencies.py
```

Control that the way the script checks for the dependency is the same name. In the default installation you'll see, for example, `bedtools`. If your installation is called `Bedtools`, or `bedtools_2.29`, or any other name, simply edit it in this script so that the tool can see it.

Also, make sure that the tool is in the `$PATH`.

#### If not in the path?

If it isn't in your `$PATH` and you don't want it to be, then simply add the full path in there rather than just the script name.

Then, add the full path also inside `nextflow.config`.


### Augustus did not get installed well

Installing augustus manually can be hard in certain systems. Did you try to install it with conda? It just takes:

```
conda install -c bioconda augustus
```

This will also install `libboost` properly, which is hard to install on your own.

#### Where is the augustus folder if I install it with Conda?

You will find the augustus main directory inside the `pkgs` subdirectory of the `conda` main directory.


## Nextflow errors

Nextflow errors usually arise from lack of internet connection of from a wrong `java` version. Make sure your connection works *from the terminal you use* and that the java version is at least 8.


## Running errors

None documented yet.


## Parameter-training errors

None documented yet.
