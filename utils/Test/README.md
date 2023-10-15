# Test
These codes provide tests of EXP-related infrastructure

### expyaml

`expyaml` reads an EXP YAML configuration file and generates a lint
report.  The code runs your configuration through the parser used by
EXP internally and will report any syntax errors.  It also checks that
mandatory map entries are present, reports whether optional map
entries are missing, and lists non-EXP entries.

Note: non-EXP entries are allowed.  For example, one can include
various bits of history and specific project info in a 'comments'
section.

### testBarrier

EXP uses a wrapper class for MPI_Barrier to help detect and report
deadlocks.  At this point, there are no known deadlocks or race
conditions.  This code demonstrates the wrapper class.
