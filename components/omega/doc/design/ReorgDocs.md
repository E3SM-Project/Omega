(omega-design-reorg-docs)=
# Documentation Reorganization

**Table of Contents**
1. [Overview](#1-overview)
2. [The starting point](#2-the-starting-point)
3. [What other ocean models do](#3-what-other-ocean-models-do)
4. [Audiences](#4-audiences)
5. [Principles](#5-principles)
6. [Design](#6-design)
7. [Scope of the initial reorganization](#7-scope-of-the-initial-reorganization)
8. [Discussion and open questions](#8-discussion-and-open-questions)
9. [Verification](#9-verification)

## 1 Overview

This document describes the reorganization of the Omega documentation from a
layout that mirrors the source tree into one organized around the questions
readers bring to it. It records the reasoning behind the new structure, the
principles that should guide future changes to the documentation, and the
decisions about mechanisms (configuration reference generation, bibliography,
versioned publication) that the structure depends on.

Like every design document in this directory, it captures our thinking at the
time it was written. It is not a description of the documentation as it exists
today; that description is the documentation itself, and the living statement of
the rules is the Developer Guide page on documentation.

## 2 The starting point

Before this reorganization, the documentation consisted of three flat lists in
the top-level `index.md`:

- a *User's Guide* of 39 pages,
- a *Developer's Guide* of 44 pages, and
- 38 *Design documents*.

The User's and Developer's Guides were organized one page per source module:
`Config`, `Broadcast`, `Halo`, `Decomp`, `Field`, `Reductions`, `TimeMgr`, and
so on, with the same list repeated in both guides. This was a deliberate crutch
while the model was being built: a developer adding a module knew exactly which
two pages to add. It has several costs that grew as the model matured.

- **Readers do not think in modules.** A user asking "how do I control the time
  step?" has to know that the answer is split across `TimeMgr`, `TimeStepping`
  and `Driver`. A developer asking "how do I add a tracer?" has to know to look
  in `Tracers`, then `Field`, then `IOStreams`.
- **The user/developer split produced duplication, not differentiation.** Many
  User's Guide pages were thinner copies of the Developer's Guide page (for
  example `Config`), and about a third of them said, in effect, "there are no
  user-configurable options; see the Developer's Guide."
- **The science had no home.** The governing equations lived only in a frozen
  design document; discretization details were scattered through
  developer pages and user pages (the second-order tracer advection derivation
  was in the *User's* Guide). Nothing in the layout told a scientist where to
  look, and the [E3SM documentation](https://docs.e3sm.org/E3SM/) pattern of a
  *Technical Guide* had no counterpart.
- **No onboarding path.** The user Quick Start was a stub, there were no
  tutorials, and "how to run Omega" was not a page.
- **No reference material.** Configuration options were documented by hand,
  inconsistently, on whichever page a developer chose; there was no single list.
- **No versioning.** The published site tracked `develop` only. With a first
  release approaching, users of a release would have had no matching
  documentation.

## 3 What other ocean models do

We surveyed the documentation of several community ocean models to find
practices worth adopting. The observations below are about *organization*, not
content.

**MOM6** ([mom6.readthedocs.io](https://mom6.readthedocs.io/en/main/)) is
organized by topic, science first: Equations, Spatial Discretization, Time
Discretization, Tracers, Grids, Parameterizations, Other Physics, Working with
MOM6, Forcing, Parallel Implementation, Testing, then an API reference and a
bibliography. Each topic has one home containing the continuous equations, the
discrete form, and pointers to code and parameters, so there is no duplication
by construction. Runtime parameters are documented automatically from the code.
The weakness is onboarding: "how do I build and run" is delegated to a wiki,
and the practical material is buried in chapter 8.

**MITgcm** ([mitgcm.readthedocs.io](https://mitgcm.readthedocs.io/en/latest/))
is a single numbered manual: Overview (theory), Discretization and Algorithm,
Getting Started, Tutorial Example Experiments, Contributing, Software
Architecture, Packages, Utilities, References. It is thorough and the tutorial
experiments are the most praised part of it. The chapters are very long.

**NEMO** ([sites.nemo-ocean.io/user-guide](https://sites.nemo-ocean.io/user-guide/))
splits cleanly into a practical *User Guide* (Basics: install, prepare, run;
Advanced: coupling, assimilation, new configurations; Appendix: release notes,
migration; plus a short Developer guide) and a separate *Reference Manual* for
the science and numerics. The user guide explicitly says it is practical rather
than theoretical and points to the reference manual for the rest.

**ROMS** (the UCLA fork on
[readthedocs](https://cworthy-ucla-roms.readthedocs.io/en/latest/); the
official documentation is a wiki) is Introduction, Key technical information,
Getting Started, Customizing, Tutorials, References, Release notes: a simple
progression from "what is it" to "how do I change it."

**Oceananigans**
([clima.github.io/OceananigansDocumentation](https://clima.github.io/OceananigansDocumentation/stable/))
is Quick start, Examples, Workflows (Models, Simulations), Concepts (Physics,
Numerical implementation), Developer, Appendix. It separates "how to use it"
from "how it works" cleanly and keeps pages short and heavily cross-linked.

**E3SM component documentation** ([docs.e3sm.org](https://docs.e3sm.org/E3SM/))
uses a *User Guide / Developer Guide / Technical Guide* triad for every
component. EAMxx is the most developed example: its Developer Guide has a
quick start, code organization, style guide, testing, "important tools and
objects" and how-tos; its Technical Guide holds physics descriptions. In other
components the Technical Guide is thin or a pointer to a PDF.

The **[Diátaxis](https://diataxis.fr)** framework is the general form of what
these sites converge on: documentation serves four distinct needs (tutorials,
how-to guides, reference, explanation), and mixing them on one page serves none
of them well.

What we take from each:

| Source | Adopted |
| ------ | ------- |
| E3SM components | The User / Developer / Technical Guide triad as the outer frame, so Omega is recognizable next to its siblings. |
| MOM6 | The science-first topical spine *inside* the Technical Guide; auto-generated parameter reference; a real bibliography. |
| NEMO | A User Guide that is explicitly practical and defers theory elsewhere. |
| MITgcm, ROMS | A place for tutorials, even if it starts empty. |
| Oceananigans | Short, focused, cross-linked pages; explicit how-to recipes for developers. |
| EAMxx | The shape of the Developer Guide: getting started, code organization, contributing, the framework objects, how-tos. |
| Diátaxis | The rule that each section answers one kind of question. |

## 4 Audiences

In priority order **at the time of writing** (before the v1 release,
planned for fall 2026):

1. **Omega and Polaris developers** adding features to the model and its test
   framework. This is the dominant audience today.
2. **Standalone users** running Omega through Polaris.
3. **E3SM coupled-model users** configuring the ocean component.
4. **E3SM integration and release staff** building and testing Omega inside
   E3SM.
5. **Scientists and reviewers** who want the equations, discretization and
   parameterizations.

This ordering will invert over time: as the model stabilizes the emphasis
shifts toward the broader E3SM development team and then the community. The
structure is chosen so that shift does not require another reorganization: the
Technical Guide exists from the start even though it begins mostly as
placeholders, and the User Guide is organized by task so that coupled-model and
community users slot in without displacing anything.

## 5 Principles

These are the rules the documentation should follow, and the rules a reviewer
should hold a documentation change to. They are restated (and maintained) in
the Developer Guide; this is the record of why they were chosen.

1. **Organize by the reader's question, not by the code.** The User Guide
   answers "how do I set up and run Omega, and what can I control?" The
   Developer Guide answers "how do I work on the code?" The Technical Guide
   answers "what does the model compute, and why?" The design documents answer
   "what were we thinking when we built it?" A page belongs where its question
   belongs, regardless of which source file implements it.
2. **One home per fact; everywhere else, link.** Every configuration option,
   equation, class and procedure is described in exactly one place. Other pages
   link to it rather than restating it. This applies across projects too: Omega
   links to Polaris and E3SM documentation rather than duplicating them, and
   vice versa.
3. **Docs ship with the code.** A pull request that changes behavior updates
   the affected guide pages in the same pull request. There is no separate
   documentation backlog.
4. **Design documents are frozen.** A design document records the reasoning at
   the time of design. It is never reconciled with the code afterward; the
   guides describe the code as it is. Design documents are labelled as such, and
   new designs still go there, using the template.
5. **Generate reference material; write everything else.** Lists that must be
   complete and exact (configuration options; eventually the API) come from a
   single machine-readable source and are generated at build time. Prose is
   written by hand.
6. **Math is written to be checked against code.** Technical Guide notation
   stays close to the variable names in the source, and each page links to the
   files that implement it, so a reviewer can verify an implementation against
   its description.
7. **Placeholders are explicit, never silent.** If a page has not been written,
   a stub page says what belongs there and links the issue that tracks it. A
   missing topic should be visible, not discovered.
8. **Every page has a stable anchor.** Anchors are prefixed by guide
   (`omega-user-`, `omega-dev-`, `omega-tech-`, `omega-design-`) so that
   cross-references survive moves within a guide and are unambiguous if the
   Omega documentation is ever merged with other E3SM documentation.
9. **Prefer short, focused pages over long chapters.** A page covers one topic
   and links to its neighbors. Readers arrive from search and from links, not
   from the beginning of a chapter.
10. **Cite the literature through a bibliography.** Technical Guide pages cite
    papers with `{cite}` roles against a single BibTeX file rather than inline
    links, so citations are consistent and collected in one place.
11. **The structure has room to grow toward science without another
    reorganization.** See [Audiences](#4-audiences).

## 6 Design

### 6.1 Top-level structure

The documentation has four top-level sections. Each answers one kind of
question and has a landing page that says who it is for and where the other
questions are answered.

| Section | Question | Primary audiences | Diátaxis form |
| ------- | -------- | ----------------- | ------------- |
| User Guide | How do I set up and run Omega, and what can I control? | 2, 3, 4 | tutorials, how-to, reference |
| Developer Guide | How do I work on the code? | 1, 4 | how-to, reference |
| Technical Guide | What does the model compute, and why? | 5, 1 | explanation |
| Design documents | What were we thinking when we built it? | 1 | (archive) |

This is a hybrid: the E3SM triad on the outside, the MOM6 topical spine on the
inside of the Technical Guide. We considered two alternatives.

*A pure audience split* (the triad with nothing more said) is what we had, and
it produced duplication rather than differentiation, because "user" and
"developer" are not different *questions* about the same topic unless the
guides are given distinct jobs. The table above gives them distinct jobs.

*A pure topical spine* (MOM6's shape) would have put the model's physics and
numerics at the top of the table of contents, with building, running,
contributing and the infrastructure classes as afterthoughts. For the next
year or two that is backwards for our audiences, and Omega's substantial
infrastructure layer (configuration, I/O streams, fields, halos,
decomposition, time management, Kokkos loops) does not fit a science spine at
all. We keep the spine, but inside the Technical Guide.

The cost of the hybrid is three trees to keep coherent. Principle 2 (one home
per fact) is what controls that cost: a configuration option's meaning lives on
one User Guide page, the discretization it selects lives on one Technical Guide
page, the class that reads it lives on one Developer Guide page, and each links
to the others.

Matching the E3SM triad exactly was judged to matter little in itself. Omega
keeps its own Sphinx build and will be linked from the E3SM site as a subpage;
the triad was adopted because it is a good fit for the questions, not for
conformity.

### 6.2 User Guide

Organized by task. Each subdirectory has a short landing page.

- **Quick start** (placeholder until the model is ready for outside users).
- **Building**: standalone build, build inside E3SM, and the CMake options
  (precision, architecture, threading, vector length, MPI on device) in one
  place.
- **Running**: what a standalone run directory needs and where Polaris takes
  over; running as the E3SM ocean component (`user_nl_omega`, compsets, how
  the coupler controls time); the outputs of a run (log files, timing files,
  error messages and what they mean).
- **Configuration**: the YAML file and how `buildnml` layers it in a coupled
  case; time management (calendar, start, stop, duration); time stepping
  choices; parallel layout (decomposition and I/O tasks); and the generated
  **configuration reference** (§6.7).
- **Input and output**: I/O streams; mesh-file requirements; initial
  conditions and restarts; available fields, groups and metadata (including
  fill values and inactive layers); analysis and diagnostics.
- **Model options**: one page per user-facing choice about the physics and
  numerics (tendency terms, tracers, advection, pressure gradient, vertical
  coordinate, vertical mixing, equation of state, forcing), each stating what
  can be set and linking the Technical Guide for what it means.
- **Tutorials** (placeholder). Standalone tutorials belong mostly to Polaris
  and will be linked; coupled E3SM tutorials may live here.

Pages of the old User's Guide whose content was "no user-configurable options"
(Broadcast, Reductions, Halo, MachEnv, Dimension, Field internals, Driver
internals, Auxiliary state, Tendencies container, Ocean state, Horizontal
operators, Tridiagonal solvers, Surface coupling internals) do not survive as
user pages; whatever user-relevant sentences they contained move to the task
page that needs them.

The User Guide owns "running Omega in E3SM" for now because no other
documentation does. It may migrate to coupled-model documentation later.

### 6.3 Developer Guide

Organized by what a developer is trying to do. Modeled on EAMxx's Developer
Guide and Oceananigans' Developer section.

- **Getting started**: quick start, conda environment, the CMake build
  system, and a new **code organization** page describing `src/`, `test/`,
  `configs/` and `cime_config/`.
- **Contributing**: linting and style, testing (unit tests and the Polaris
  `omega_pr` suite), documentation (this section restates the principles and
  says where each kind of content goes, how to build the docs, and how to write
  a design document), and CIME integration (`buildnml`).
- **Framework**: the infrastructure classes a developer uses, grouped by
  concern on the landing page: parallel programming (data types, Kokkos loops,
  machine environment, decomposition, halos, broadcasts, reductions); runtime
  services (configuration, logging, errors, timers, time manager); fields and
  I/O (dimensions, fields and metadata, parallel I/O, streams); numerical
  utilities (tridiagonal solvers).
- **Model**: how the ocean model is put together: the drivers, state,
  mesh and operators, vertical coordinate, auxiliary variables, tendencies,
  time steppers, tracers, equation of state, pressure gradient, vertical
  advection and mixing, forcing, surface coupling, analysis. These pages
  describe classes and their use, not the mathematics; each links to its
  Technical Guide counterpart.
- **How-to recipes**: short task pages ("add a configuration option", "add a
  tracer", "add a field to output", "add a tendency term", "add an analysis
  operator", "add a time stepper", "add a supported mesh", "add a unit test").
  This is the Diátaxis how-to form and the thing developers *using* the guide
  most often want. A few are written from existing material; the rest are
  placeholders.

The Developer Guide no longer carries the old rule that it must describe the
mathematical terms and their discretized form. That content moves to the
Technical Guide (§6.4). A reviewer verifying an implementation reads the
Technical Guide page and the Developer Guide page; we accepted that cost for
the sake of a single home for the math.

The Framework and Model pages still correspond roughly to classes, which is
appropriate for reference material about classes. What has changed is that they
are grouped by concern, they no longer have a user-guide twin, and they no
longer carry the science.

### 6.4 Technical Guide

MOM6's spine, adapted to Omega, and thin at first:

- **Equations**: notation; governing equations; the vertical coordinate
  (pseudo-height).
- **Discretization**: the MPAS mesh specification; TRiSK horizontal operators
  and discrete auxiliary quantities (vorticity, kinetic energy); tracer
  advection (horizontal and vertical); pressure gradient; time stepping
  (including split-explicit); tridiagonal systems.
- **Physics and parameterizations**: equation of state; vertical mixing;
  surface forcing; parameterizations.
- **Diagnostics**: definitions of analysis quantities.
- **Bibliography**.

Existing mathematical content is moved here from the old user and developer
pages (the mesh specification, second-order tracer advection, pressure-gradient
background, vertical-coordinate definition, vertical mixing formulas, equation
of state details, tridiagonal system forms, kinetic energy on cells). New
science content — in particular a proper governing-equations chapter derived
from the V1 design document — is placeholders with tracking issues.

Technical Guide pages follow principle 6: variable names close to the code, a
link to the implementing files, and `{cite}` references.

### 6.5 Design documents

Unchanged in content. They get a landing page stating that they are frozen
records of design-time reasoning (principle 4) and that the guides describe the
code as it is. The template stays and is referenced from the Developer Guide's
documentation page.

### 6.6 Directory layout, naming and anchors

```
doc/
  index.md              landing: what Omega is, how the docs are organized
  userGuide/            index.md + build/ run/ config/ io/ model/ tutorials/
  devGuide/             index.md + gettingStarted/ contributing/ framework/ model/ howto/
  techGuide/            index.md + equations/ discretization/ physics/ + Bibliography.md, references.bib
  design/               index.md + existing documents, Template.md
```

Each subdirectory has an `index.md` landing page with a toctree; the sidebar
shows two levels. File names remain CamelCase for consistency with the
existing files. Anchors on every page use the prefixes from principle 8; the
old `omega-` (no guide) prefix on some user pages is normalized to
`omega-user-`.

Old URLs are not preserved. At the time of the reorganization nothing outside
the repository links to individual pages except Polaris, which is updated in
tandem, and maintaining redirects would have been ongoing cost for no reader.

### 6.7 Configuration reference

The configuration reference is generated (principle 5) from two files that
parallel each other:

- `configs/Default.yml`, which already exists, holds the default *value* of
  every option and is what `omega_buildnml` validates a case's configuration
  against.
- `configs/ConfigDescriptions.yml`, a new companion file with the same nesting
  as `Default.yml`, holds the *meaning* of every option: a description, and
  optionally a type, units and the allowed values or range. Each section also
  carries a description and the anchor of the User Guide page that discusses
  it.

For example, the `Decomp` section of `Default.yml`,

```yaml
Decomp:
  HaloWidth: 3
  DecompMethod: MetisKWay
```

would be paired with

```yaml
Decomp:
  description: Horizontal domain decomposition across MPI tasks.
  doc: omega-user-config-parallel
  HaloWidth:
    description: Number of halo layers around each subdomain. Must be at least 3 for higher-order tracer advection.
    type: int
    range: [1, null]
  DecompMethod:
    description: Partitioning algorithm.
    type: string
    values: [MetisKWay, ParMetisKWay]
```

A Python script, run by Sphinx at the start of every documentation build,
reads both files and writes the reference page: one section per top-level
configuration group, one table per section with the option name, default,
type, units, allowed values and description. The generated page is not
committed; it is rebuilt every time.

The script also checks the two files against each other, and any failure is a
build error, so the pull-request documentation build fails (principle 7):

- every option in `Default.yml` has a description;
- every entry in `ConfigDescriptions.yml` names an option that exists in
  `Default.yml`;
- every default value satisfies its own `type`, `values` or `range`;
- every `doc` anchor exists.

This is what keeps the reference from rotting: a pull request that adds an
option to `Default.yml` cannot pass CI without describing it, and one that
removes or renames an option cannot leave a stale description behind.

The descriptions file is also useful beyond the documentation build, because
it is ordinary YAML that any tool can read. `omega_buildnml` can use `values`
and `range` to reject a `user_nl_omega` override that names a real option but
gives it an illegal value, at `case.setup` time rather than at run time. The
model itself could do the same check at startup later.

Three approaches were considered:

1. *Descriptions as comments in `Default.yml`*: no new file, but comments have
   no structure to check against and are easy to leave stale.
2. *A companion descriptions file* (chosen): explicit and checkable, at the
   cost of a second file that the build check keeps in step with the first.
3. *MOM6-style descriptions in the C++ `Config::get` calls*, dumped by a
   documentation mode of the executable: keeps each description next to the
   code that reads the option, but requires touching every module and building
   the executable to build the documentation.

Option 3 remains attractive and is an open discussion item (§8); the generator
is written so that its input could later come from such a dump instead of the
companion file.

Only the design, the placeholder page and the tracking issue are part of the
initial reorganization; the descriptions file, generator and build check
follow.

### 6.8 Bibliography

`sphinxcontrib-bibtex` is added to the documentation build with a single
`techGuide/references.bib`. Technical Guide pages cite with `{cite:t}` and
`{cite:p}`; `techGuide/Bibliography.md` renders the list. Design documents are
frozen and keep whatever citation style they have.

### 6.9 Versioned publication

The published site keeps one directory per version, following the Polaris
pattern:

- The publish workflow runs on pushes to `develop`, on published releases, and
  on manual dispatch. It builds into `_build/html/<version>` where `<version>`
  is the branch name or tag (`DOCS_VERSION`), replaces only that directory on
  `gh-pages` under `Omega/`, never touches other versions, and regenerates
  `Omega/shared/versions.json` from the directories present.
- A small JavaScript version switcher in the sidebar reads `versions.json`.
  `develop` is listed first, then release versions sorted numerically.
- `Omega/index.html` redirects to `develop/`.

Rebuilding old versions' documentation from within a newer checkout was
rejected from experience: it is fragile and rarely worth it. Old versions are
simply kept as built.

### 6.10 Placeholders and tracking issues

A placeholder page consists of the page's anchor and title, a standard
admonition saying the page has not yet been written, a sentence or two on what
it will contain, and a link to the tracking issue. Placeholders are listed in
the tracking issues so they can be found and retired.

### 6.11 Coordination with Polaris

Standalone-run documentation is shared between the two projects by linking:
Polaris documents the workflow (setting up tasks, building, running suites);
Omega documents what is Omega-specific (configuration, streams, mesh and
initial-condition requirements, outputs). A Polaris pull request accompanies
the reorganization to update any links into the Omega documentation and to add
links to the new pages.

## 7 Scope of the initial reorganization

In scope for the reorganization pull request:

- the new directory structure, landing pages and toctrees;
- moving and retitling every existing page into its new home, merging the
  user/developer duplicates and splitting mathematical content into the
  Technical Guide;
- the new Developer Guide pages that are connective (code organization,
  documentation rules) and the how-to recipes that can be assembled from
  existing text;
- the design documents landing page;
- the bibliography mechanism;
- versioned publication;
- explicit placeholders, with tracking issues, for everything below.

Deferred to tracked issues:

- a real user Quick Start, running standalone and running in E3SM;
- the configuration descriptions file, generator and build check;
- migrating the governing equations and other science into the Technical
  Guide;
- the remaining how-to recipes;
- tutorials;
- splitting the longest existing pages (principle 9).

## 8 Discussion and open questions

- **Descriptions in code (option 3 of §6.7).** Keeping the description next
  to the `Config::get` call that reads an option is the most robust way to keep
  the two in step, at the cost of a documentation mode in the executable. Worth
  revisiting once the companion-file approach is in place and its maintenance
  cost is known.
- **Where coupled-run documentation ultimately lives.** Owned by Omega for now
  (§6.2); a coupled-model documentation space in E3SM may be the right home
  later.
- **Tutorials: here or Polaris.** A Tutorials section exists in the User Guide
  as a placeholder. Standalone tutorials will most likely be Polaris pages that
  Omega links; coupled tutorials may be written here.
- **API reference.** MOM6 and Oceananigans generate one. Omega has no Doxygen
  setup; whether an API reference earns its place is left open.

## 9 Verification

- The documentation builds with `make html-strict` (warnings are errors), as
  enforced on every pull request.
- Every page has an anchor with the correct prefix; no page from the old layout
  is left orphaned or unlinked.
- Every placeholder links a tracking issue.
- The versioned publication workflow is exercised by a `develop` build and, at
  the first release, by a tag build, with the switcher listing both.
- The Polaris documentation's links into Omega resolve.
