# TOC graphic

From the package root:

```bash
python -B 03_figures/toc/make_figure.py --output ../toc_result/TOC_graphic.png
```

This redraws the accepted TOC schematic from Python drawing objects. It does
not copy the final manuscript PNG. No trajectories, analysis CSVs, ChimeraX or
files outside this folder are needed as inputs. The output must be a new file
outside the package. Only one PNG is exported, at 975 × 525 pixels and 300 dpi.

## Files

- `make_figure.py`: safe public entry, dependency/font checks and final-image check.
- `make_toc_graphic_schematic.py`: original drawing functions.
- `make_toc_graphic_schematic_v1_refined.py`: accepted layout refinements and pixel checks.
- `TOC_graphic_schematic.png`: baseline image used only to validate the original drawing.
- `TOC_graphic.png`: accepted final reference image, also displayed in the repository README.
- `source_manifest.json`: source hashes, exact environment and accepted final PNG hash.

The two drawing modules are copies of the manuscript project's TOC sources.
Their drawing functions are unchanged. Direct execution is disabled in these
copies so that old output paths and the legacy PNG/TIFF exporter cannot be used
accidentally. Run `make_figure.py` instead.

The tested environment is recorded in `source_manifest.json`. DejaVu Sans is
selected explicitly to match the accepted artwork, even if other system fonts
are installed. Any baseline-pixel or final-file mismatch stops the run; the
entry does not silently change the design or substitute another image.

The final PNG included here as `TOC_graphic.png` is an
independent copy with the same expected hash. It is not read by this renderer.
`03_figures/run.py` continues to generate the 15 numbered figures and six tables;
the TOC is a separate, lightweight entry and not a sixteenth numbered figure.

## Historical review presentations

[review/](review/README.md) retains the Chinese and English seven-slide PPTX
presentations from the 2026-09-04 design review. Both contain their own images
and can be opened independently of the original design directories. They are
historical comparison material, not the final approved TOC or drawing inputs.
The accepted V1 refinement remains the output of `make_figure.py` above.

中文：用上面的命令单独重绘摘要图。这里包含完整绘图依赖，基准图只用于校验，
不是最终图的替代来源。程序不会覆盖所交付的参考图片，也不会写入原始实验目录。
