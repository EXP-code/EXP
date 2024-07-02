# Paper

This is the JOSS abstract/paper describing the EXP release

## Contents

All of the text is in the `paper` directory.  Contents are:

| File      | Description |
| ---       | ---         |
| paper.md  | The primary markdown text file |
| paper.bib | The bibliography file |
| paper.pdf | The compiled PDF file |


## Compile to PDF with Docker as follows

```
docker run --rm --volume $PWD/paper:/data --user $(id -u):$(id -g) --env JOURNAL=joss openjournals/inara
```
