---
title: Methodology
---

The content presented in this blog is compiled from deep dives into practical underlying workflows, combined with my interpretations of the algorithmic and biological concepts described in published research.

## Seurat architecture and execution

Seurat v5.5.0 is structured as an R-centric workflow where the core logic and primary user interfaces are written in R. However, to handle heavy computational tasks efficiently, Seurat and its dependencies leverage C, C++, and Fortran routines. C and Fortran functions can be called natively within R, whereas C++ functions are integrated via the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) interface.

## Conceptual Dissection

To fully comprehend an implementation, I dissect it into two complementary components:
* **The "Body":** What specific execution steps occur internally?
* **The "Spirit":** Why is each step designed and performed this way?

I use interactive debuggers to explore the "body" and map statistical and mathematical concepts from articles, documentation, or lectures to capture the "spirit". When a paper's method section is particularly thorough, it significantly accelerates the process of dissecting the codebase.

## Debugging Workflow

To inspect the functional commands in the tutorial, I use the native R debugger from the CLI:

| Command | Purpose |
| :--- | :--- |
| `<R> debug(<func>)` | Set a persistent breakpoint on a function |
| `<R> undebug(<func>)` | Remove a persistent breakpoint from a function |
| `<R> browser()` | Insert a manual hard-coded breakpoint into source code |
| `<R> n` (or Enter) | Execute next line (step over) |
| `<R> s` | Step into a function call |
| `<R> Sys.getpid()` | Get current R process ID (PID) for `gdb` attachment |
| `<R> saveRDS(<var list>, file = <path>)` | Save current frame variables to an RDS file |

The corresponding R source code acts as a visual guide during execution manually. When execution transitions to compiled external code (C/C++/Fortran), I use `gdb` to attach to the R process by PID and step through the source code similarly.

::: {note}

Target shared (C/C++/Fortran) libraries must be built in debug mode from source. You can then run a new R script reproducing the target command by loading the variable stack and the debuggable shared library. Example R scripts are available in the [commands](https://github.com/Truongphi20/seurat_blog/tree/main/commands) repository folder.

:::

| Command | Purpose |
| :--- | :--- |
| `<bash> cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build` | Build source binaries with debug symbols |
| `<bash> gdb -p <PID>` | Attach GDB to running R process |
| `<gdb> break <pos>` | Set breakpoint by function name or file position |
| `<gdb> c` | Continue process execution |
| `<gdb> s` | Step into a function |
| `<gdb> dir <path>` | Map source code directory paths |

In cases where logic resides within core R packages (such as `stats`), building `r-base` from source with debug symbols enabled is necessary (see [setup-r.sh](https://github.com/Truongphi20/seurat_blog/blob/main/.devcontainer/setup-r.sh)).

Key call-stack breakpoints for each functional step are mapped out in [breakpoint diagrams](https://github.com/Truongphi20/seurat_blog/blob/main/docs/static/breakpoints).


## Environment

A [DevContainer](https://code.visualstudio.com/docs/devcontainers/containers) based on the `satijalab/seurat:5.5.0` base image (contains R v4.5.2) serves as the reproducible development environment for code exploration and site deployment. Detailed setup configs can be found in the [Dockerfile](https://github.com/Truongphi20/seurat_blog/blob/main/.devcontainer/Dockerfile).