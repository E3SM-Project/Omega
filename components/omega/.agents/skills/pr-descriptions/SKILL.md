---
name: pr-descriptions
description: Write or update the description of an Omega pull request. Use when opening a pull request on E3SM-Project/Omega or editing its body.
---

# Pull request descriptions

The reader is deciding whether to review. The description is also the
message of the merge commit, which ends up in E3SM's `master`, so it
follows E3SM's
[Commit and PR description template](https://e3sm.atlassian.net/wiki/spaces/DOC/pages/16253965/Commit+and+PR+description+template)
and the pull request template.

- The pull request title is the commit title: imperative, under 70
  characters, the whole change, no issue numbers.
- The opening part is one short paragraph of a few sentences: what
  changed and why. Not how. Start with a verb ("Add", "Fix"), not "This
  PR". Everything else goes below `---`.
- Say how users are affected: a changed interface, config option or
  answers.
- End the opening part with `[BFB]`, `[non-BFB]`, `[NCC]` or `[CC]`, plus
  `[NML]` if a config option is added, removed or changes its default.
  Apply the matching labels.
- Wrap the part above `---` at about 72 characters, and use Markdown
  there as anywhere else (backticks for file names, for example). Do not
  hard-wrap anything below `---`.
- Below a line holding only `---`: `Fixes` lines and other issue
  references, links, figures and tables. Name the repository in every
  reference (`E3SM-Project/Omega#123`).
- Anything needing a reviewer decision goes in its own short list below
  `---`, never mid-paragraph.
- A list of changed behaviours is fine. A trace of the mechanism is not.
- No commit list. No testing; that goes in a separate `Testing` comment.
- Several fixes usually means several pull requests.

## Calibration

Measured over Omega and Polaris pull requests merged in 2023 and 2024,
before any agent wrote there. Descriptions run 36 words at the median, 66
at the seventy-fifth percentile, 114 at the ninetieth, and 398 at the
longest. In E3SM over the same years, the part above `---` runs 54 words
at the median and 142 at the ninetieth percentile. A recent agent-written
Polaris description ran 1167 words.

## Enough

A real Omega description (E3SM-Project/Omega#131), as it reads under the
rules above:

> The `updateTimeLevels` routine copied the whole state between the host
> and the device twice. Reuse `exchangeHalo` instead, which fixes this
> performance issue, and remove some unused variables.
>
> [BFB]

## Too much

A real agent-written description gave each of five fixes its own section
and traced its mechanism:

> `polaris.ocean.conservation` gained a module-scope import of
> `polaris.ocean.model.time`, which runs
> `polaris/ocean/model/__init__.py`, which imports the step classes, which
> import `polaris.ocean.conservation` back. The module could no longer be
> imported on its own; it worked only when something else imported
> `polaris.ocean.model` first.

Someone deciding whether to review does not need the cycle traced. One
sentence would do; the rest belongs in the commit message.
