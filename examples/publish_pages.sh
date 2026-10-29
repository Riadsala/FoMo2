#!/bin/sh
##########################################################################
# Publish the rendered example pages to GitHub Pages.
#
# The .html pages are not tracked on main (each re-render would add
# another ~30 MB to the history). Instead this script replaces the
# gh-pages branch with a single commit holding the current pages, plus an
# index page listing them, so only the latest version is ever kept.
#
# Render the pages you changed first, then, from anywhere in the repo:
#   sh examples/publish_pages.sh
# The site is https://riadsala.github.io/FoMo2/ (GitHub: Settings ->
# Pages -> deploy from branch gh-pages, once).
##########################################################################

set -eu

root=$(git rev-parse --show-toplevel)
remote=$(git -C "$root" remote get-url origin)
site=$(mktemp -d)
trap 'rm -rf "$site"' EXIT

cd "$root"
pages=$(find examples -name '*.html' -not -path '*_files/*' | sort)
[ -n "$pages" ] || { echo "No rendered pages in examples/: render them first." >&2; exit 1; }

# section <level> <heading> <folder>: a heading (<h2>, <h3>, ...), then
# each of the folder's pages, copied to the same path on the site and
# listed under it
section() {
  echo "<h$1>$2</h$1>"
  echo "<ul>"
  for page in $(echo "$pages" | grep "^$3/"); do
    mkdir -p "$site/$(dirname "$page")"
    cp "$page" "$site/$page"
    title=$(sed -n 's:.*<title>\(.*\)</title>.*:\1:p' "$page" | head -1)
    echo "<li><a href=\"$page\">${title:-$page}</a> <code>${page#$3/}</code></li>"
  done
  echo "</ul>"
}

{
  cat <<'EOF'
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>FoMo 2.0 examples</title>
<style>
  body { font-family: system-ui, sans-serif; max-width: 42rem; margin: 2rem auto; padding: 0 1rem; line-height: 1.5; }
  code { font-size: 0.9em; color: #555; }
</style>
</head>
<body>
<h1>FoMo 2.0: worked examples</h1>
<p>Rendered pages from the <a href="https://github.com/Riadsala/FoMo2">FoMo 2.0 repository</a>.
Start with the single-level spatial example; <a href="https://github.com/Riadsala/FoMo2/blob/main/examples/README.md">examples/README.md</a> gives the reading order.</p>
EOF
  echo "<h2>Simulation examples</h2>"
  section 3 "Single-level models" examples/simulated/single_level
  section 3 "Multi-level models" examples/simulated/multi_level
  section 3 "Model comparison" examples/simulated/model_comparison
  section 2 "Fitting FoMo to empirical datasets" examples/real_data
  cat <<'EOF'
</body>
</html>
EOF
} > "$site/index.html"
touch "$site/.nojekyll"  # serve the files as they are, without Jekyll

cd "$site"
git init -q -b gh-pages
git add -A
git commit -q -m "Publish example pages"
git push -f "$remote" gh-pages
echo "Published $(echo "$pages" | wc -l | tr -d ' ') pages to gh-pages."
