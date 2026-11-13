# Omega Agent Instructions

These instructions apply to everything under `components/omega` and to
pull requests on E3SM-Project/Omega. They add to the root `AGENTS.md`,
and take precedence over it where the two differ.

## Source of truth

- Style is enforced by `pre-commit`, configured in the root
  `.pre-commit-config.yaml`, with `.clang-format`, `ruff.toml`, `mypy.cfg`
  and `.fortitude.toml` in this directory. If anything here conflicts with
  these tools, follow the tools.
- The Developer's and User's Guides in `doc/` (published at
  <https://docs.e3sm.org/Omega/>) describe how Omega works. Read the
  relevant page before changing a part of the code, and link to it rather
  than explaining unchanged features again.

## Environment

- Linting and the docs build use the conda environment described in
  `doc/devGuide/CondaEnv.md`. Run `pre-commit` on the files you changed
  before finishing, and fix problems rather than suppressing them.
- Omega is built and tested on supported machines. Ask the developer
  before submitting batch jobs, which spend a shared allocation.

## Building and testing

- What a pull request needs is set by the pull request template and
  `doc/devGuide/Testing.md`: CTests, the Polaris `omega_pr` suite and
  coupled CIME tests, on at least one CPU and one GPU machine.
- Polaris has optional helpers. Use them only when the developer wants
  them:
  - `polaris setup` and `polaris suite` build Omega themselves; see the
    [Polaris quick start](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html).
  - `utils/omega/ctest` builds Omega and runs the CTests. See the
    `omega-ctest` skill.
  - `utils/omega/pr_testing` is one way to run everything the template
    asks for. See the `omega-pr-testing` skill.
  - `utils/bisect` finds the first commit at which a Polaris task fails or
    changes answers.
- Never tick checklist boxes in a pull request on the developer's behalf.

## Documentation

- A new feature updates the User's Guide and Developer's Guide, and adds
  its pages to `doc/index.md`.
- Start a design document from `doc/design/Template.md`.

## Branches

- Name branches `omega/<feature>`, on your own fork of E3SM. The form
  `<user>/omega/<feature>` is allowed but discouraged, since the fork
  already says whose branch it is.
- Pull requests target `develop`.

## Commit messages and pull request descriptions

Omega's history, merge commits included, becomes part of E3SM's `master`.
A pull request's description becomes the message of its merge commit, so
both follow E3SM's
[Commit and PR description template](https://e3sm.atlassian.net/wiki/spaces/DOC/pages/16253965/Commit+and+PR+description+template).
In short:

- The title is imperative ("Fix", "Add"), under 70 characters, and
  describes the whole change. It never contains an issue or pull request
  number.
- There is always a body, even if it only repeats the title. A few
  sentences say what changed and why, for a developer reading `git log`
  years from now who has not seen the pull request. Say it in the code's
  own terms, and say how users are affected, such as a changed interface
  or changed answers.
- Wrap commit messages, and the part of a pull request description above
  `---`, at about 72 characters. Do not hard-wrap anything below `---`.
- No URLs in the body, other than DOIs.

Omega adds:

- Name the repository in every reference to an issue or pull request:
  `E3SM-Project/Omega#123`, `E3SM-Project/E3SM#9071`,
  `E3SM-Project/polaris#897`. Once a commit reaches E3SM, a bare `#123`
  links to E3SM's own #123. Closing keywords accept the long form
  (`Fixes E3SM-Project/Omega#123`). Commit hashes do not change on the
  merge into E3SM, so they need no qualifier.
- Do not use issue, pull request, work-package or design-decision numbers
  as names for things. A reference may sit beside a description, but the
  message must make sense without following it.

A pull request description also:

- Opens with one short paragraph of a few sentences. Everything else
  goes below `---`.
- Ends its opening part with one of `[BFB]`, `[non-BFB]`, `[NCC]` or
  `[CC]`, plus `[NML]` if it adds or removes a config option or changes a
  default. The pull request gets the matching labels.
- Puts issue references, `Fixes` lines, links, figures and tables below a
  line holding only `---`.
- May use Markdown in its opening part. This overrides the root
  `AGENTS.md`, which asks for plain text there.
- Keeps only the template checklist lines that apply, and lists neither
  commits nor testing. Testing goes in a separate `Testing` comment.

## GitHub comments and issues

- Do not hard-wrap comments or issues. GitHub wraps them for display, and
  hard breaks make later edits show up as reflowed paragraphs.
- An issue says what happened, what was expected instead, and enough
  about the configuration and commands to reproduce it.

## Writing for human readers

These rules apply to anything a colleague reads: pull request
descriptions, GitHub comments, issues, plans. Not code comments or commit
messages, where a reader who wants the mechanism is already in the right
place. Per-artifact rules and worked examples are in
`.agents/skills/<artifact>/SKILL.md` in this directory, as plain Markdown.
Agents that do not load skills by themselves should read the matching one
before writing. (`.claude/skills` is a link to the same directory.)

Assume the reader is an Omega developer who knows the overall structure
and purpose of the code. Explain what changed, particularly conceptual
changes, and nuance that informs the discussion. Do not explain unchanged
features again.

Write less; do not pack the same content into denser sentences. Keep
headings, tables and links. The problem is length.

- **Lead with the answer.** The first two sentences say what you found,
  changed or propose. Setup and reproduction go last.
- **One point per paragraph, and few paragraphs.** Colleagues write one
  to three per comment. Say each thing once.
- **Do not narrate the mechanism.** The chain of calls, and why the fix
  is right, go in the commit message. Here, say what broke and where to
  look.
- **Cut clauses that qualify rather than inform**, and any sentence whose
  only job is to justify the one before it.
- **Use backticks about half as often as feels natural.** They are for
  what a reader would type or grep.
- **One document, one decision.** Anything still relevant after a pull
  request merges belongs in an issue, not a comment.

Sign anything posted on the user's behalf: pull request descriptions,
issues and comments. Put the signature at the end of a pull request
description, below the `---` line, and on its own line at the top of
anything else:

```
*Posted by <agent> on @<user>'s behalf*
```

If the text grew out of the user's own draft or substantial direction,
use instead:

```
*Drafted by @<user> and <agent> together*
```

Name the agent, not the vendor: `Claude Code`, `Codex`, and so on.
