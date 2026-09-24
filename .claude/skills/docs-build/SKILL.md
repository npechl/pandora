---
name: docs-build
description: Build and sanity-check Pandora's documentation site (Zensical, config at zensical.toml, source under docs/) and apply this project's doc-writing conventions. Use this whenever you finish editing any file under docs/, whenever the user asks to "build the docs", "check the docs site", "does the site still build", or mentions broken links/anchors in the documentation, and before telling the user a docs change is done — a docs edit isn't finished until it's been built.
---

# Building Pandora's docs

The site is built with Zensical. The config is `zensical.toml` at the repo root, source pages are under `docs/`, and the output goes to `site/`.

To draft or edit a new page from scratch with the user, use `doc-coauthoring`, then come back here to apply any conventions. Use this skill for refining and checking existing pages.

## 1. Build

From the repo root:

```bash
uv run zensical build --clean --strict
```

`--strict` makes the build exit non-zero on any warning (a broken link or a missing anchor). Without it the build exits 0 even with warnings. A non-zero exit means the change isn't done: fix every warning it prints. If `zensical` isn't found, run `uv sync --all-extras` first.

## 2. Checks Zensical doesn't do

Zensical passes both of these problems without a warning, so run both checks after every build.

**Duplicate anchors.** Two headings with the same explicit `{: #id }` produce duplicate HTML ids. Deep links then land on the first one.

```bash
for f in $(find site -name '*.html'); do grep -o ' id="[^"_][^"]*"' "$f" | sort | uniq -d | sed "s|^|$f:|"; done
```

Any output is a clash to fix. The `[^"_]` skips the `__codelineno-*` ids that mkdocstrings repeats on purpose.

**Pages missing from `nav`.** Navigation is an explicit list in `zensical.toml`. A page left out of it builds fine but can't be reached from the site.

```bash
uv run python -c "import pathlib,re; nav=set(re.findall(r'\"([^\"]+\.md)\"', pathlib.Path('zensical.toml').read_text())); print([str(p.relative_to('docs')) for p in pathlib.Path('docs').rglob('*.md') if str(p.relative_to('docs')) not in nav])"
```

It should print `[]`. Otherwise, add each listed page to `nav`.

## 3. Conventions

**One line per paragraph.** Write each paragraph and each list or glossary entry as a single line, however long. Don't hard-wrap prose at 80 characters; that limit is for Python code only. Leave code blocks and tables as they are. When you edit a paragraph on a page that is still hard-wrapped, join the paragraph you touched into one line.

**Reference pages: headings and stable anchors.** Group terms under `##` sections. Give each term its own `###` heading with an explicit anchor, so the URL survives renames and other pages can deep-link to it:

```markdown
### Term name {: #term-slug }
```

When a `##` section and a `###` term inside it have the same text, give the section its own id, e.g. `## Canonicalisation {: #canonicalisation-section }`. Otherwise the second heading is silently renamed `…_1`.

**Cross-link.** Link related terms to each other by their anchors. Link from reference content to the matching `docs/usage/*.md` guide, `docs/recipes/*.md` recipe or `docs/reference/*.md` page, so a reader always has somewhere to go next.
