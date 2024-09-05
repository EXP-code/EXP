# User modules

All sub directories will be scheduled for building by CMake.

This allows you to include your own code and especially modular force
routines into EXP simulations. We suggest following the examples
`src/user` to get started.  The following steps will get you started
quickly:

1. Make a sub directory of your own choosing and set that as your
   working directory
2. Copy the `CMakeLists.txt` file from `src/user/CMakeLists.txt`
3. Pick one of the modules for your template.  For example, if you are
   interested a applying an external force, check out the
   `UserMNdisk.cc` and header `UserMNDisk.H`.  This produces the force
   and potential Miyamoto-Nagai disk.  Copy these to your new
   sub directory.  Rename them to something mnemonic and adjust they to
   suit your needs.
4. Edit `CMakeLists.txt` to include the name of your new module and
   source code.  Don't forget to remove the original libraries from
   `src/user` that you are not using.
5. Then, your next compile will compile and install your code as where
   it can be found and used by EXP.

You can add any code to the directory and it will be built and linked
to EXP libraries.  This is useful for building standalone applications
that use EXP classes or pyEXP for your own intricate analysis. For
this, use the code in `utils` as a pattern for your own applications.

