# MicroclimateTalk

Slides and demo code for a 10-minute talk on [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) and [BiophysicalGrids.jl](https://github.com/BiophysicalEcology/BiophysicalGrids.jl), given at Montpellier 2026.

Built with Quarto + reveal.js. Figures are regenerated from Julia in CI and the deck is published to GitHub Pages.

## Local build

```bash
julia --project --threads=auto src/generate_figures.jl
quarto render
```

Open `_site/talk.html`.

## Layout

- `talk.qmd` — slide deck
- `src/point_demo.jl` — single-point microclimate (Montpellier)
- `src/grid_demo.jl` — gridded microclimate (Mont Aigoual, Cévennes)
- `src/generate_figures.jl` — driver that runs both demos
- `figures/` — generated PNG/GIF outputs (not committed)
- `.github/workflows/build.yml` — render + deploy to GitHub Pages
