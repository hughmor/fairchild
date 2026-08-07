# Logos

| File | Ink | Use on |
|---|---|---|
| `logo.svg` | `#20242b` | **light** backgrounds |
| `logo_dark.svg` | `#edebe4` | **dark** backgrounds |
| `logo_icon.svg` | `#20242b` on a `#edebe4` panel | square mark, light |
| `logo_icon_dark.svg` | `#edebe4` on a `#20242b` panel | square mark, dark |

The names say what the file is *for*, not what colour it is — `logo_dark.svg` is
the light-coloured one, because it goes on a dark page. Easy to get backwards;
if the wordmark ever looks washed out, that is what happened.

The full logos are transparent. The icons carry their own background panel, so
they are the ones to use where a mark needs to sit on an unknown colour: GitHub
repo avatar, social preview card, a favicon.

## Palette

| | |
|---|---|
| Ink | `#20242b` |
| Cream | `#edebe4` |
| Amber | `#e9a63c` (`#f0b352` on dark) |

## Switching automatically

Markdown that has to work in both GitHub themes needs a `<picture>`, because a
plain `![](…)` cannot respond to the theme:

```html
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/logos/logo_dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="docs/logos/logo.svg">
  <img alt="fairchild" src="docs/logos/logo.svg" width="440">
</picture>
```

The `<img>` fallback should be the **light-background** file: it is what renders
where `prefers-color-scheme` is unavailable, which is most non-browser Markdown
viewers, and those are overwhelmingly light.
