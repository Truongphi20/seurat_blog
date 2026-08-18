---
title: Methodology
---

The content presented in this blog is compiled from deep-dives into the practical workflows underlying, alongside my interpretations of the algorithmic and biological concepts found in published researches.

## Seurat architecture and execution

Seurat v5.5.0 works as an R-centric workflow, where whole dependencies and internal logics are written in R. Moreover, Seurat and its dependences employ C/C++/Fortran-function calls to tackle heavy computational tasks. Specificly, C or Fortran functions are able to be instinctively called in R code directly while C++ functions go through declaration of the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) tunnel.

## Reverse-engineering

To fully comprehence an implementation, I dissect it into two parts: the "body" - What are implemented steps inside?; and the "spirit" - Why does each step being performed like this?

I utilize programmatic debugger to explore the "body", and map the statictical/mathematical concepts from the articles, documentation, or lectures to completely capture the "spirit".

### Identifying functional steps by debugging 

In general, I use R debuger directly on command line to walk into functional commands in the tutorial, used commands:

| Command           |   Purpose     |
|:-------------     | :----------   |
| `<R> debug(<func>)`   |  Set persistent breakpoint on a function      |
| `<R> undebug(<func>)`	|  Remove persistent breakpoint from a function |
| `<R> browser()`	    |  Insert into code as a manual hard-coded breakpoint |
| `<R> n` (or Enter)	|  Execute next line (step over) |
| `<R> s`	            |  Step into a function call     |
| `<R> Sys.getpid()`    |  Get current R PID on computer (to get hijacked by `gdb`) | 
| `<R> saveRDS(<var list>, file = <save path>)` | Store temporary variables on current stack to RDS file |             

Noticeably, corresponded R source code is seemed as the manually line tracker while in debuging process. When external functions (C/C++/Fortran) are called, `gdb` is employed to hijack the process by its PID and walking through the source code simmilarly. 

Attentionally, used binary shared library should be build in debug mode from its source code beforehand (use `cmake`). Then, running a new R script, which captures the specific investigated command, with loading frame and debuggable shared library. See R scripts in my [commands](https://github.com/Truongphi20/seurat_blog/tree/main/commands) folder for illustration.

| Command           |   Purpose     |
|:-------------     | :----------         |
| `<bash> cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build`  | Build binary source code with debug mode  |
| `<bash> gdb -p <PID>`    | Hijacking R process |
| `<gdb> break <pos>`     | Set breakpoints by function name or position in source code |
| `<gdb> c`               | Continue the process |
| `<gdb> s`               | Step into a function |
| `<gdb> dir`             | Mapping path of source code folder | 

Occasionally, the methodology hides behind the R foundational packages, e.g. `stats`, the build of r-base with debuggable mode is neccessary (see [setup-r.sh](https://github.com/Truongphi20/seurat_blog/blob/main/.devcontainer/setup-r.sh)). By default, debuggable r is setup when igniting the [devcontainer](https://code.visualstudio.com/docs/devcontainers/containers).

Important breakpoints for functional steps are diagramed as stack flow charts (see [breakpoint diagrams](https://github.com/Truongphi20/seurat_blog/blob/main/docs/static/breakpoints)).

### Conceptual Dissection

## Environment

A DevContainer based on the `satijalab/seurat:5.5.0` image is used as the environment (obtains R v4.5.2) for walking through the codebase and deploying this blog. Detailed configurations are available in the [Dockerfile](https://github.com/Truongphi20/seurat_blog/blob/main/.devcontainer/Dockerfile).