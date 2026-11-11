---
name: review-comments
description: Write a review comment, review findings, or a reply to review feedback on an Omega pull request. Use when reviewing code, reporting what testing someone else's branch turned up, or answering a reviewer's question.
---

# Review comments

The reader is deciding what to change.

- Put each finding as an inline comment on the line it concerns, one
  point each. That is where colleagues put them, and it is why their
  review bodies are short.
- The review body summarizes: what you ran, and the verdict. Two or three
  sentences.
- Use a list in the body only for requests that span files.
- No section on what already works. One line for all of it, if any.
- Say what you could not check.

## Calibration

Measured over Omega and Polaris pull requests from 2023 and 2024, before
any agent wrote there, and leaving out pasted output. Review bodies run
21 words at the median and 69 at the ninetieth percentile. Inline
comments run 16 at the median, 57 at the ninetieth percentile and 329 at
the longest. E3SM's are a little shorter. A recent agent-written review
ran 1117 words, with the findings starting 444 words in.

## Enough

Inline, one point and a suggestion:

> I think this else branch is not needed. I believe that
> `create_mirror_view_and_copy` already does the right thing (does not
> copy) when the execution space passed to it is accessible.

> Intel compiler (correctly) flagged this as a duplicate of the same
> function in Decomp, so we can either rename them, protect them with a
> namespace or break these out into a utility module somewhere.

A review body that says what was run and the verdict:

> Verified that unit tests pass and that halo tests work with different
> numbers of mpi tasks on pm-cpu. The changes look good.

In the body, when several requests span the whole change:

> Thanks for putting this together and the overall system looks good.
> However, I think we need to make this a cleaner, simpler PR with the
> following changes:
>
> - For the CIME changes and YAKL submodule, I think we need to sync the
>   OMEGA repo so it's up to date with E3SM.
> - Can you please remove the logger and spdlog. This will need a separate
>   PR and discussion.

## Too much

A real agent-written review spent its first 444 words on "How this was
reviewed", "What the previous review asked for" and four paragraphs of
"What works", then traced each finding's mechanism:

> It is static: computed once in the `MOC` constructor from `NumBins`,
> `MinLat` and `MaxLat`, and never updated. But it is attached to the
> output streams with `addField()` like any other field, so every
> reduction in every file carries a copy of the same 61 numbers.
