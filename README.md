# CollectiveCellMigration
Analysis of collective cell migration, including PIV, cellular orientation analysis, cell division analysis, etc

Version 3 - 21-10-2022
Changes:
  - Changed the way how cell divisions are detected. It now uses a combination of segmentation via uNet, followed by peakfinding in the "cross correlation space" and subsequent classification. Changes result in both increased accuracy/robustness and faster analysis

Version 2 - 29-08-2022
Changes:
- Added sample images together with results, plots, etc.
- Improved speed
- Added some checks for correct parameter input
- Improved internal readability
- Added additional .pdf as documentation

Version 1 - 30-03-2022
Features:
- Automatic image reading with subsequent PIV, orientation and cell division analysis for phase contrast images of dense cell layers
- Analysis tools for velocity and orientation fields
- Plotting of output data
