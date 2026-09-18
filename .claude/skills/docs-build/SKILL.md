---
name: docs-build
description: Build and sanity-check Pandora's documentation site (Zensical, config at zensical.toml, source under docs/) and apply this project's doc-writing conventions. Use this whenever you finish editing any file under docs/, whenever the user asks to "build the docs", "help me build the docs for a specific page", "check the docs site", "does the site still build", or mentions broken links/anchors in the documentation, and before telling the user a docs change is done — a docs edit isn't finished until it's been built.
---

# Building Pandora's docs

Pandora's documentation site is built with Zensical (a Material-for-MkDocs-style static site generator). Config lives in `zensical.toml` at the repo root; source pages live under `docs/`.

## 1. Build and read the output

From the repo root:

```bash
zensical build
```

Treat any warning or error in the output as something to fix before calling the docs change done — this is the fastest way to catch broken internal links, duplicate heading anchors, and other issues that are easy to introduce by hand and easy to miss by eye. A clean build prints something like `No issues found`.

If `zensical` isn't on `PATH`, say so rather than skipping the check silently — don't guess that the docs are fine.

## 2. Conventions to apply while editing

These conventions were settled on for this project; apply them without being asked so the docs stay consistent across pages.

**Don't hard-wrap prose.** Markdown treats a single `\n` inside a paragraph as a soft break, so hard-wrapping lines at ~80 chars (the Python line-length convention from `pyproject.toml`) buys nothing visually and only makes docs harder to diff and grep. Write each paragraph or list/glossary entry as one line, however long.

**Group reference-style content under headings, and give terms stable anchors.** For glossary-like or reference pages, group related terms under `##` section headings rather than one flat list — `zensical.toml` already enables `toc` with `permalink = true`, so headings automatically populate the page's in-page table of contents, giving readers a real way to scan and jump around instead of just scrolling. Give each individual term its own `###` heading with an explicit anchor via the `attr_list` extension (also already enabled):

```markdown
### Term name {: #term-slug }
```

An explicit anchor keeps the URL stable even if the heading text changes later, and lets other pages deep-link straight to that term.

Watch for anchor collisions: if a `##` section heading and a `###` term inside it share the same text (e.g. a "Canonicalisation" section containing a "Canonicalisation" term), they'll auto-slug to the same id. Give the section heading its own explicit id in that case (e.g. `## Canonicalisation {: #canonicalisation-section }`) — `zensical build` will flag the clash as a duplicate-anchor warning if you miss it.

**Cross-link to encourage exploration.** Link related terms to each other using their anchors, and link out from reference content to the matching `docs/usage/*.md` guide, `docs/recipes/*.md` recipe, or `docs/reference/*.md` API page. The goal is that a reader following one link keeps finding another one worth clicking, rather than dead-ending on a single page.

**Update `nav` in `zensical.toml` for new pages.** The site's navigation is an explicit list in `zensical.toml`, not derived from the `docs/` directory structure. A new page that isn't added to `nav` will build without error but simply won't show up anywhere in the site — check for this whenever a docs edit adds a new file.
